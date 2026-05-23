import numpy as np
import scipy.io
import scipy.linalg
import sklearn.metrics
from torch import nn


def kernel(ker, X1, X2, gamma):
    K = None
    if not ker or ker == 'primal':
        K = X1
    elif ker == 'linear':
        if X2 is not None:
            K = sklearn.metrics.pairwise.linear_kernel(
                np.asarray(X1).T, np.asarray(X2).T)
        else:
            K = sklearn.metrics.pairwise.linear_kernel(np.asarray(X1).T)
    elif ker == 'rbf':
        if X2 is not None:
            K = sklearn.metrics.pairwise.rbf_kernel(
                np.asarray(X1).T, np.asarray(X2).T, gamma)
        else:
            K = sklearn.metrics.pairwise.rbf_kernel(
                np.asarray(X1).T, None, gamma)
    return K


class TCA:
    def __init__(self, kernel_type='primal', dim=30, lamb=1, gamma=1):
        '''
        Init func
        :param kernel_type: kernel, values: 'primal' | 'linear' | 'rbf'
        :param dim: dimension after transfer
        :param lamb: lambda value in equation
        :param gamma: kernel bandwidth for rbf kernel
        '''
        self.kernel_type = kernel_type
        self.dim = dim
        self.lamb = lamb
        self.gamma = gamma

    def fit(self, Xs, Xt):
        '''
        Transform Xs and Xt
        :param Xs: ns * n_feature, source feature
        :param Xt: nt * n_feature, target feature
        :return: Xs_new and Xt_new after TCA
        '''
        X = np.hstack((Xs.T, Xt.T))
        # Safe column norm — cells with zero expression across all genes
        # must not produce NaN (a pre-existing bug exposed by full data).
        norms = np.linalg.norm(X, axis=0)
        X[:, norms > 0] /= norms[norms > 0]
        m, n = X.shape
        ns, nt = len(Xs), len(Xt)

        if self.kernel_type == 'primal':
            # ── Memory-efficient TCA ──────────────────────────────────
            # Original formulation constructs two (n×n) matrices M and H,
            # which cost O(n²) memory (~25 GB for 56K samples).
            #
            # M = e·eᵀ / ‖e·eᵀ‖_F  (rank-1, where e = [1/ns; -1/nt])
            # H = I - (1/n)·1·1ᵀ    (centering matrix)
            #
            # For primal kernel K = X  →  a, b are (m, m) = (541, 541):
            #
            #   K·M·Kᵀ = (X·e)·(X·e)ᵀ / (1/ns+1/nt)    [vector outer product]
            #   K·H·Kᵀ = X·Xᵀ - (1/n)·(X·1)·(X·1)ᵀ     [gram − rank-1 update]
            #
            # Memory: O(m·n + m²) ≈ 244 MB instead of O(n²) ≈ 50 GB.
            # ──────────────────────────────────────────────────────────

            # X·e  — difference of column means (m, 1)
            Xe = X[:, :ns].sum(axis=1, keepdims=True) / ns \
               - X[:, ns:].sum(axis=1, keepdims=True) / nt
            norm_factor = 1 / ns + 1 / nt            # ‖e·eᵀ‖_F
            a = (Xe @ Xe.T) / norm_factor + self.lamb * np.eye(m)

            # X·Xᵀ and X·1  for centering (m, m) and (m, 1)
            XXt = X @ X.T
            X1 = X.sum(axis=1, keepdims=True)
            b = XXt - (X1 @ X1.T) / n

            K = X  # primal: kernel is identity
        else:
            # Fallback — for linear / rbf kernels K is (n × n) so the
            # (n × n) M and H matrices are unavoidable here.
            K = kernel(self.kernel_type, X, None, gamma=self.gamma)
            e = np.vstack((1 / ns * np.ones((ns, 1)), -1 / nt * np.ones((nt, 1))))
            M = (e * e.T) / (1 / ns + 1 / nt)
            H = np.eye(n) - np.ones((n, n)) / n
            n_eye = n
            a = K @ M @ K.T + self.lamb * np.eye(n_eye)
            b = K @ H @ K.T

        w, V = scipy.linalg.eig(a, b)
        ind = np.argsort(w)
        A = V[:, ind[:self.dim]]
        Z = A.T @ K
        norms_z = np.linalg.norm(Z, axis=0)
        Z[:, norms_z > 0] /= norms_z[norms_z > 0]

        Xs_new, Xt_new = Z[:, :ns].T, Z[:, ns:].T
        return Xs_new, Xt_new


activation_dict = {
    'relu': nn.ReLU,
    'sigmoid': nn.Sigmoid
}


class sample_MLPEncoder(nn.Module):
    """sample MLPEncoder encoder model for scSpace."""

    def __init__(self, input_size, common_size, hidden_size, activation):
        """Init MLP encoder."""
        super(sample_MLPEncoder, self).__init__()
        self.restored = False
        common_size = common_size
        self.linear = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            activation_dict[activation](),
            nn.Linear(hidden_size, common_size)
        )

    def forward(self, x):
        out = self.linear(x)
        return out


