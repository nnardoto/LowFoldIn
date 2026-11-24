"""PBE (sem spin) em malha periódica: V_xc = ∂e/∂ρ − 2∇·(∂e/∂σ ∇ρ), σ = |∇ρ|².
Derivadas parciais de e(ρ, σ) por passo complexo (exatas à precisão de máquina); gradiente/divergência por FFT.
Também PW92 puro (LDA) para comparação. Unidades: Hartree, bohr."""
import numpy as np
KAPPA, MU = 0.804, 0.2195149727645171
BETA = 0.06672455060314922; GAMMA = (1 - np.log(2)) / np.pi ** 2
def _pw92(rs):
    A, a1, b1, b2, b3, b4 = 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / (2 * A * (b1 * np.sqrt(rs) + b2 * rs + b3 * rs ** 1.5 + b4 * rs ** 2)))
def e_xc(rho, sig, gga=True):
    """densidade de energia xc por volume (Ha/bohr³)."""
    kF = (3 * np.pi ** 2 * rho) ** (1 / 3)
    ex = -3 * kF / (4 * np.pi)
    rs = (3 / (4 * np.pi * rho)) ** (1 / 3)
    ec = _pw92(rs)
    if not gga:
        return rho * (ex + ec)
    s2 = sig / (2 * kF * rho) ** 2
    Fx = 1 + KAPPA - KAPPA / (1 + MU * s2 / KAPPA)
    ks = np.sqrt(4 * kF / np.pi)
    t2 = sig / (2 * ks * rho) ** 2
    Ac = BETA / GAMMA / (np.exp(-ec / GAMMA) - 1)
    At2 = Ac * t2
    Hc = GAMMA * np.log(1 + BETA / GAMMA * t2 * (1 + At2) / (1 + At2 + At2 ** 2))
    return rho * (ex * Fx + ec + Hc)
class Grid:
    def __init__(self, A, n):
        self.A = np.asarray(A, float); self.n = tuple(n)
        Bm = 2 * np.pi * np.linalg.inv(self.A).T
        m = [np.fft.fftfreq(n[k], 1.0 / n[k]) for k in range(3)]
        M = np.stack(np.meshgrid(*m, indexing="ij"), -1)
        # zera a frequência de Nyquist (derivada ímpar mal definida)
        nyq = np.zeros(n, bool)
        for k in range(3):
            if n[k] % 2 == 0: nyq |= (M[..., k] == -n[k] // 2)
        self.G = M @ Bm; self.G[nyq] = 0.0
    def grad(self, f):
        F = np.fft.fftn(f)
        return np.stack([np.fft.ifftn(1j * self.G[..., k] * F).real for k in range(3)], -1)
    def div(self, v):
        return sum(np.fft.ifftn(1j * self.G[..., k] * np.fft.fftn(v[..., k])).real for k in range(3))
def vxc(rho, grid, gga=True, rmin=1e-10):
    rho = np.maximum(rho, rmin)
    g = grid.grad(rho); sig = np.sum(g ** 2, -1)
    h = 1e-30
    de_dr = e_xc(rho + 1j * h * rho, sig, gga).imag / (h * rho)
    if not gga:
        return de_dr
    hs = h * np.maximum(sig, 1e-30)
    de_ds = e_xc(rho.astype(complex), sig + 1j * hs, gga).imag / hs
    return de_dr - 2 * grid.div(de_ds[..., None] * g)
class GridFD(Grid):
    """diferenças centrais de 2ª ordem ao longo dos eixos da malha (como a malha do OpenMX), periódicas."""
    def __init__(self, A, n):
        self.A = np.asarray(A, float); self.n = tuple(n)
        self.Jinv = np.linalg.inv(self.A / np.array(n)[:, None])   # ∇f = Jinv @ (∂f/∂i)
    def _di(self, f):
        return np.stack([(np.roll(f, -1, k) - np.roll(f, 1, k)) / 2 for k in range(3)], -1)
    def grad(self, f):
        return self._di(f) @ self.Jinv.T
    def div(self, v):
        # v = Σ_k c_k ∂r/∂i_k = J^T c  →  c = J^{-T} v (componentes contravariantes); ∇·v = Σ_k ∂c_k/∂i_k (volume da célula constante)
        c = v @ self.Jinv
        return sum((np.roll(c[..., k], -1, k) - np.roll(c[..., k], 1, k)) / 2 for k in range(3))
