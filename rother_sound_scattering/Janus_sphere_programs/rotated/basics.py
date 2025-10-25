#            Modul "basics.py"
#
# This module containes:
# 
# subroutines to calculate the:  
#
# separation matrix "S(l,n,nue,x)" needed to transform
# the regular and outgoing wave functions if the coordinate 
# system is shifted along the z-axis
#
# the matrix of rotation to transform the regular and outgoing wave 
# functions if the coordinate system is rotated by use of the 
# Eulerian angles of rotation
#
# the C21- and Q12 operators of the block operator for arbitrarily  
# oriented bispheres (for both the rotation-shift, and the translation 
# formulation)
#
# the T-matrices of acoustically soft (Dirichlet problem), hard 
# (von Neumann problem), and penetrable (mixed problem) spheres
#
# the matrices "Y_l_J_u" and "Y_l_J_l" related to the Janus spheres.  
# The elements are calculated from integrals over products of Legendre 
# polynomials
#
# the diagonal matrices required for the Janus sphere
#
# the T-matrix for the s-p-Janus sphere
#
# the T-matrix for the h-p-Janus sphere
#
# the T-matrix for the h-s-Janus sphere
#
# the T-matrix for the sp-2layer sphere
#
# the T-matrix for the hp-2layer sphere
#


import numpy as np
import scipy.misc as scm
import scipy.special as scs
import scipy.integrate as sqa


# Calculation of the separation matrix "S(l,n,nue,x)" according to 
# [Martin, S. 104, Gl. (3.126)]:

def SM(el, en, nu, x):
    zi = 0 + 1j
    n1 = np.absolute(el) 
    w = zi / 2 / x
    we = np.exp(zi * x) / zi / x
    pf = (-1)**n1 * zi**(en + nu) * we
    j_u = n1
    j_o = en + nu + 1
    S_sumj = 0.
    for j in range(j_u, j_o):
        du = [0, j - en]
        do = [nu, j - n1]
        s_u = max(du)
        s_o = min(do) + 1
        j1 = scs.factorial(j, exact=True)
        j2 = w**j
        S_sums = 0.        
        for s in range(s_u, s_o):
            jms = j - s
            f1 = np.sqrt(2. * nu + 1)
            f1 = np.log(f1)
            f2 = np.sqrt(2. * en + 1)
            f2 = np.log(f2)
            f = f1 + f2
            a1 = 0.
            for i in range(2, nu + s + 1):
                a1 = a1 + np.log(float(i))
            a2 = 0.
            for i in range(2, en + jms + 1):
                a2 = a2 + np.log(float(i))
            a = a1 + a2
            b1 = 0.
            for i in range(2, s + 1):
                b1 = b1 + np.log(float(i))
            b1s = 0.
            for i in range(2, n1 + s + 1):
                b1s = b1s + np.log(float(i))
            b2 = 0.
            for i in range(2, nu - s + 1):
                b2 = b2 + np.log(float(i))
            b3 = 0.
            for i in range(2, jms + 1):
                b3 = b3 + np.log(float(i))
            b3s = 0.
            for i in range(2, jms - n1 + 1):
                b3s = b3s + np.log(float(i))
            b4 = 0.
            for i in range(2, en - jms + 1):
                b4 = b4 + np.log(float(i))
            b = b1s + b1 + b2 + b3s + b3 + b4
            c1 = 0.
            for i in range(2, nu + n1 + 1):
                c1 = c1 + np.log(float(i))                    
            c2 = 0.
            for i in range(2, en - n1 + 1):
                c2 = c2 + np.log(float(i))
            c = c1 + c2
            d1 = 0.
            for i in range(2, nu - n1 + 1):
                d1 = d1 + np.log(float(i))
            d2 = 0.
            for i in range(2, en + n1 + 1):
                d2 = d2 + np.log(float(i))
            d = d1 + d2
            s_fak = np.exp(c - d)
            s_fak = np.sqrt(s_fak)
            s_pref = np.log(s_fak)
            res_fak = np.exp(a + f + s_pref - b)
            S_sums = S_sums + res_fak
        S_sumj = S_sumj + j1 * j2 * S_sums
    S_M = pf * S_sumj
    return S_M    


# Calculation of the matrix of rotation "drm(n,l,l1,alpha,theta_p,gamma)"
# (note, that "alpha", "theta_p", and "gamma" must be given in degree!):

def drm(en, el, el_1, alpha, theta_p, gamma):
    alpha = alpha * np.pi / 180.
    gamma = gamma * np.pi / 180.
    zi = 0 + 1j
    n_p_l1 = en + el_1
    n_m_l1 = en - el_1
    n_p_l = en + el
    n_m_l = en - el
    l_p_l1 = el + el_1
    s0 = [0, - l_p_l1]
    s1 = [n_m_l, n_m_l1]
    s_min = max(s0)
    s_max = min(s1)
    z1 = 0.
    for i in range(2, n_p_l1 + 1):
        z1 = z1 + np.log(float(i))
    z1 = z1 * 0.5
    z2 = 0.
    for i in range(2, n_m_l1 + 1):
        z2 = z2 + np.log(float(i))
    z2 = z2 * 0.5
    z3 = 0.
    for i in range(2, n_p_l + 1):
        z3 = z3 + np.log(float(i))
    z3 = z3 * 0.5
    z4 = 0.
    for i in range(2, n_m_l + 1):
        z4 = z4 + np.log(float(i))
    z4 = z4 * 0.5
    za = z1 + z2 + z3 + z4
    D_m = 0.
    for s in range(s_min, s_max + 1):
        n1 = 0.
        for i in range(2, s + 1):
            n1 = n1 + np.log(float(i))
        n2 = 0.
        for i in range(2, l_p_l1 + s + 1):
            n2 = n2 + np.log(float(i))
        n3 = 0.
        for i in range(2, n_m_l1 - s + 1):
            n3 = n3 + np.log(float(i))
        n4 = 0.
        for i in range(2, n_m_l - s + 1):
            n4 = n4 + np.log(float(i))
        na = n1 + n2 + n3 + n4
        d_m = np.exp(za - na)
        d_m = d_m * (-1.0)**(n_p_l1 + s)
        e_p0 = 2 * en - 2 * s - l_p_l1
        e_p1 = 2 * s + l_p_l1
        w_sin = (scs.sindg(theta_p / 2.)**e_p0 if e_p0 != 0 else 1.0)
        w_cos = scs.cosdg(theta_p / 2.)**e_p1
        d_m = d_m * w_sin * w_cos
        D_m = D_m + d_m
    D_m = np.exp(- zi * el * alpha) * D_m * np.exp(- zi * el_1 * gamma)
    return D_m


# Calculation of the C21- and Q12 operators of the block operator for 
# arbitrarily oriented bispheres!

def Q12(n1cut, alpha, theta_p, sm_l_i_j, t_o, qr):
    c_l0_m0 = []
    for l0 in range(-n1cut, n1cut + 1):    
        c_m0_3 = []
        for m0 in range(np.abs(l0), n1cut + 1):
            c1 = t_o[m0]
            c2 = 0.0
            for l1 in range(-n1cut, n1cut + 1):
                q_l1 = qr[l1 + n1cut]
                z_l1 = sm_l_i_j[np.abs(l1)]
                z_l1_m0 = z_l1[m0]
                D_rueck = drm(m0, l0, l1, alpha, theta_p, 0.)
                for n1 in range(np.abs(l1), n1cut + 1):
                    n1z = n1 - np.abs(l1)
                    S_rueck = z_l1_m0[n1]
                    c2 = c2 + (-1)**m0 * q_l1[n1z] * S_rueck * D_rueck
            c2 = c2 * c1 * (-1)
            c_m0_3 = c_m0_3 + [c2]
        c_l0_m0 = c_l0_m0 + [c_m0_3]
    return c_l0_m0    


def Q12_tl_it1(n1cut, alpha, theta_p, sm_l_i_j, t_o, qr):
    c_l0_m0 = []
    for l0 in range(-n1cut, n1cut + 1):    
        c_m0_3 = []
        for m0 in range(np.abs(l0), n1cut + 1):
            c1 = t_o[m0]
            c2 = 0.0
            for l1 in range(-n1cut, n1cut + 1):
                q_l1 = qr[l1 + n1cut]
                z_l1 = sm_l_i_j[np.abs(l1)]
                z_l1_m0 = z_l1[m0]
                D_rueck = drm(m0, l0, l1, alpha, theta_p, 0.)
                for n1 in range(np.abs(l1), n1cut + 1):
                    D_hin = drm(n1, l1, 0, 0., -theta_p, -alpha)
                    D_ges = D_hin * D_rueck
                    n1z = n1 - np.abs(l1)
                    S_rueck = z_l1_m0[n1]
                    c2 = c2 + (-1)**m0 * q_l1[n1z] * S_rueck * D_ges
            c2 = c2 * c1 * (-1)
            c_m0_3 = c_m0_3 + [c2]
        c_l0_m0 = c_l0_m0 + [c_m0_3]
    return c_l0_m0    

    
def C21(n1cut, alpha, theta_p, sm_l_i_j, t_o1, cr):
    q_l1_n1 = []
    for l1 in range(-n1cut, n1cut + 1):
        x_l1 = sm_l_i_j[np.abs(l1)]    
        q_n1_3 = []
        for n1 in range(np.abs(l1), n1cut + 1):
            q1 = t_o1[n1]
            x_l1_n1 = x_l1[n1]
            q2 = 0.0
            for l0 in range(-n1cut, n1cut + 1):
                c_l0 = cr[l0 + n1cut]
                for m0 in range(np.abs(l0), n1cut + 1):
                    m0z = m0 - np.abs(l0)
                    S_rueck = x_l1_n1[m0]
                    D_rueck = drm(m0, l0, l1, alpha, theta_p, 0.)
                    D_rueck = np.conjugate(D_rueck)
                    q2 = q2 + (-1)**m0 * c_l0[m0z] * S_rueck * D_rueck
            q2 = q2 * q1 * (-1)
            q_n1_3 = q_n1_3 + [q2]
        q_l1_n1 = q_l1_n1 + [q_n1_3]     
    return q_l1_n1


def C21_tl_it1(n1cut, alpha, theta_p, sm_l_i_j, t_o1, cr):
    q_l1_n1 = []
    for l1 in range(-n1cut, n1cut + 1):
        q_n1_3 = []
        for n1 in range(np.abs(l1), n1cut + 1):
            q1 = t_o1[n1]
            q2 = 0.0
            for l0 in range(-n1cut, n1cut + 1):
                c_l0 = cr[l0 + n1cut]
                x_l0 = sm_l_i_j[np.abs(l0)]
                x_l0_m0 = x_l0[n1]
                D_hin = drm(n1, l1, l0, alpha, theta_p, 0.)
                for m0 in range(np.abs(l0), n1cut + 1):
                    m0z = m0 - np.abs(l0)
                    S_rueck = x_l0_m0[m0]
                    D_rueck = drm(m0, l0, 0, 0., -theta_p, -alpha)
                    D_ges = D_hin * D_rueck
                    q2 = q2 + (-1)**m0 * c_l0[m0z] * S_rueck * D_ges
            q2 = q2 * q1 * (-1)
            q_n1_3 = q_n1_3 + [q2]
        q_l1_n1 = q_l1_n1 + [q_n1_3]     
    return q_l1_n1

# Calculation of the T-matrix of an acoustically soft sphere:

def tm_s(ncut, beta):
    zi = 0. + 1.0j
    n = np.linspace(0,ncut,ncut+1,dtype=int)
    bes_1 = scs.spherical_jn(n,beta,derivative=False)
    neu_1 = scs.spherical_yn(n,beta,derivative=False)
    han_1 = bes_1 + zi * neu_1
    tm_soft = - bes_1 / han_1
    return tm_soft

    
# Calculation of the T-matrix of an acoustically hard sphere:

def tm_h(ncut, beta):
    zi = 0. + 1.0j
    n = np.linspace(0,ncut,ncut+1,dtype=int)
    bes_1 = scs.spherical_jn(n,beta,derivative=True)
    neu_1 = scs.spherical_yn(n,beta,derivative=True)
    han_1 = bes_1 + zi * neu_1
    tm_hard = - bes_1 / han_1
    return tm_hard


# Calculation of the T-matrix of an acoustically penetrable sphere:

def tm_p(ncut, beta, a, k, k_p_k, rho_p):
    zi = 0. + 1.0j
    beta_p = k * k_p_k * a
    kap_k = k_p_k
    kap_r = rho_p
    n = np.linspace(0,ncut,ncut+1,dtype=int)
    bes = scs.spherical_jn(n,beta)
    bes_1 = scs.spherical_jn(n,beta,derivative=True)
    bes_p = scs.spherical_jn(n,beta_p)
    bes_p_1 = scs.spherical_jn(n,beta_p,derivative=True)        
    neu = scs.spherical_yn(n,beta)
    neu_1 = scs.spherical_yn(n,beta,derivative=True)
    han = bes + zi * neu
    han_1 = bes_1 + zi * neu_1
    tz = kap_k / kap_r * bes * bes_p_1 - bes_1 * bes_p
    tn = kap_r / kap_k * han_1 * bes_p - bes_p_1 * han
    tm_pene = kap_r / kap_k * tz / tn
    return tm_pene
    

# Calculation of the matrices of Janus spheres:

def a_le_p(l, n, x):
    alp_ln = scs.lpmn(n, n, np.cos(x))
    alep = alp_ln[0]
    alep = alep[l]
    return alep
        
def Y_l_J_u(l, ncut, theta_j):
    al = np.abs(l)
    theta_j = theta_j * np.pi / 180.
    Y_l_J_n_m = []
    for n in range(al, ncut + 1):
        Y_l_J_m = []
        for m in range(al, ncut + 1):
            mf = scm.factorial((m - al), exact=True) / \
                    scm.factorial((m + al), exact=True)
            nf = scm.factorial((n - al), exact=True) / \
                    scm.factorial((n + al), exact=True)
            c_n_m = np.sqrt((2. *n + 1) * (2. * m + 1) * mf * nf) / 2.
            y_n_m = lambda x: np.sin(x) * a_le_p(al, ncut, x)[n] * \
                        a_le_p(al, ncut, x)[m]
            Y_l_J_m = Y_l_J_m + [c_n_m * sqa.quad(y_n_m, theta_j, np.pi, \
                                        epsabs=1.0e-04, epsrel=1.0e-04)[0]]
        Y_l_J_n_m = Y_l_J_n_m + [Y_l_J_m]
    P_l_Ju_n_m = np.array(Y_l_J_n_m)
    return P_l_Ju_n_m

def Y_l_J_l(l, ncut, theta_j):
    al = np.abs(l)
    theta_j = theta_j * np.pi / 180.
    Y_l_J_n_m = []
    for n in range(al, ncut + 1):
        Y_l_J_m = []
        for m in range(al, ncut + 1):
            mf = scm.factorial((m - al), exact=True) / \
                    scm.factorial((m + al), exact=True)
            nf = scm.factorial((n - al), exact=True) / \
                    scm.factorial((n + al), exact=True)
            c_n_m = np.sqrt((2. *n + 1) * (2. * m + 1) * mf * nf) / 2.
            y_n_m = lambda x: np.sin(x) * a_le_p(al, ncut, x)[n] * \
                        a_le_p(al, ncut, x)[m]
            Y_l_J_m = Y_l_J_m + [c_n_m * sqa.quad(y_n_m, 0.0, theta_j, \
                                        epsabs=1.0e-04, epsrel=1.0e-04)[0]]
        Y_l_J_n_m = Y_l_J_n_m + [Y_l_J_m]
    P_l_Jl_n_m = np.array(Y_l_J_n_m)
    return P_l_Jl_n_m
    
def diag_h(l, beta, ncut):
    al = np.abs(l)
    zi = 0. + 1.0j
    D_h = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta)
        neu = scs.spherical_yn(i,beta)
        han = bes + zi * neu        
        D = []
        for j in range(al, ncut + 1):
            diag = (han if j == i else 0.0)
            D = D + [diag]
        D_h = D_h + [D]
    D_h = np.matrix(D_h)
    return D_h

def diag_j(l, beta, ncut):
    al = np.abs(l)
    D_j = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta)
        D = []
        for j in range(al, ncut + 1):
            diag = (bes if j == i else 0.0)
            D = D + [diag]
        D_j = D_j + [D]
    D_j = np.matrix(D_j)
    return D_j

def diag_jp(l, beta_p, ncut):
    al = np.abs(l)
    D_jp = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta_p)
        D = []
        for j in range(al, ncut + 1):
            diag = (bes if j == i else 0.0)
            D = D + [diag]
        D_jp = D_jp + [D]
    D_jp = np.matrix(D_jp)
    return D_jp

def diag_hs(l, beta, ncut):
    al = np.abs(l)
    zi = 0. + 1.0j
    D_hs = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta,derivative=True)
        neu = scs.spherical_yn(i,beta,derivative=True)
        han = bes + zi * neu        
        D = []
        for j in range(al, ncut + 1):
            diag = (han if j == i else 0.0)
            D = D + [diag]
        D_hs = D_hs + [D]
    D_hs = np.matrix(D_hs)
    return D_hs

def diag_js(l, beta, ncut):
    al = np.abs(l)
    D_js = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta,derivative=True)
        D = []
        for j in range(al, ncut + 1):
            diag = (bes if j == i else 0.0)
            D = D + [diag]
        D_js = D_js + [D]
    D_js = np.matrix(D_js)
    return D_js

def diag_jps(l, k_p_k, rho_p, beta_p, ncut):
    al = np.abs(l)
    kappa = rho_p / k_p_k  
    D_jps = []
    for i in range(al, ncut + 1):
        bes = scs.spherical_jn(i,beta_p,derivative=True) / kappa
        D = []
        for j in range(al, ncut + 1):
            diag = (bes if j == i else 0.0)
            D = D + [diag]
        D_jps = D_jps + [D]
    D_jps = np.matrix(D_jps)
    return D_jps

def e_matrix(l, ncut):
    al = np.abs(l)
    c11 = []
    for i in range(al, ncut + 1):
        c1 = []
        for j in range(al, ncut + 1):
            c = (1.0 if j == i else 0.0)
            c1 = c1 + [c]
        c11 = c11 + [c1]
    E = np.matrix(c11)
    return E


# Calculation of the T-matrix for the h-s-Janus sphere:
    
def t_hs(l, beta, theta_j, ncut):
    m_h = diag_h(l, beta, ncut)
    m_hs = diag_hs(l, beta, ncut)
    m_j = diag_j(l, beta, ncut)
    m_js = diag_js(l, beta, ncut)
    pa_l_j = Y_l_J_l(l, ncut, theta_j)
    pb_l_j = Y_l_J_u(l, ncut, theta_j)
    
    m_h_p = pb_l_j * m_h
    m_hs_p = pa_l_j * m_hs
    m_j_p = pb_l_j * m_j
    m_js_p = pa_l_j * m_js
    A = m_hs_p + m_h_p
    B = m_js_p + m_j_p
    tm_hs = - np.linalg.inv(A) * B
    return tm_hs
        

# Calculation of the T-matrix for the s-p-Janus sphere:
    
def t_sp(l, beta, beta_p, theta_j, k_p_k, rho_p, ncut):
    m_h = diag_h(l, beta, ncut)
    m_j = diag_j(l, beta, ncut)
    m_hs = diag_hs(l, beta, ncut)
    m_js = diag_js(l, beta, ncut)
    pu_l_j = Y_l_J_u(l, ncut, theta_j)
    pl_l_j = Y_l_J_l(l, ncut, theta_j)
    d_jp = diag_jp(l, beta_p, ncut)   
    d_jps = diag_jps(l, k_p_k, rho_p, beta_p, ncut)

    m_dcs = m_hs - k_p_k / rho_p * m_h * d_jps * np.linalg.inv(d_jp)
    m_dds = m_js - k_p_k / rho_p * m_j * d_jps * np.linalg.inv(d_jp)
    A = pl_l_j * m_h + pu_l_j * m_dcs
    B = pl_l_j * m_j + pu_l_j * m_dds
    tm_sp = - np.linalg.inv(A) * B
    return tm_sp


# Calculation of the T-matrix for the h-p-Janus sphere:
    
def t_hp(l, beta, beta_p, theta_j, k_p_k, rho_p, ncut):
    m_h = diag_h(l, beta, ncut)
    m_j = diag_j(l, beta, ncut)
    m_hs = diag_hs(l, beta, ncut)
    m_js = diag_js(l, beta, ncut)
    pu_l_j = Y_l_J_u(l, ncut, theta_j)
    pl_l_j = Y_l_J_l(l, ncut, theta_j)
    d_jp = diag_jp(l, beta_p, ncut)   
    d_jps = diag_jps(l, k_p_k, rho_p, beta_p, ncut)

    m_dcs = m_h - rho_p * m_hs * d_jp * np.linalg.inv(d_jps) / k_p_k
    m_dds = m_j - rho_p * m_js * d_jp * np.linalg.inv(d_jps) / k_p_k
    A = pl_l_j * m_hs + pu_l_j * m_dcs
    B = pl_l_j * m_js + pu_l_j * m_dds
    tm_hp = - np.linalg.inv(A) * B
    return tm_hp


# Calculation of the T-matrix of the sp-2-layer sphere:

def tm_sp(ncut, beta_b, a, b, k, k_p_k, rho_p):
    zi = 0. + 1.0j
    beta_pa = k * k_p_k * a
    ts = tm_s(ncut, beta_pa)
    beta_pb = k * k_p_k * b
    kap_k = k_p_k
    kap_r = rho_p
    n = np.linspace(0,ncut,ncut+1,dtype=int)
    bes = scs.spherical_jn(n,beta_b)
    bes_1 = scs.spherical_jn(n,beta_b,derivative=True)
    bes_b = scs.spherical_jn(n,beta_pb)
    bes_1_b = scs.spherical_jn(n,beta_pb,derivative=True)       
    neu = scs.spherical_yn(n,beta_b)
    neu_1 = scs.spherical_yn(n,beta_b,derivative=True)
    han = bes + zi * neu
    han_1 = bes_1 + zi * neu_1
    neu_b = scs.spherical_yn(n,beta_pb)
    neu_1_b = scs.spherical_yn(n,beta_pb,derivative=True)
    han_b = bes_b + zi * neu_b
    han_1_b = bes_1_b + zi * neu_1_b
    bes_p = scs.spherical_jn(n,beta_pb) + ts[n] * han_b            
    bes_p_1 = scs.spherical_jn(n,beta_pb,derivative=True) + ts[n] * han_1_b
    tz = kap_k / kap_r * bes * bes_p_1 - bes_1 * bes_p
    tn = kap_r / kap_k * han_1 * bes_p - bes_p_1 * han
    tm_2l_sp = kap_r / kap_k * tz / tn
    return tm_2l_sp
    

# Calculation of the T-matrix of the hp-2-layer sphere:

def tm_hp(ncut, beta_b, a, b, k, k_p_k, rho_p):
    zi = 0. + 1.0j
    beta_pa = k * k_p_k * a
    th = tm_h(ncut, beta_pa)
    beta_pb = k * k_p_k * b
    kap_k = k_p_k
    kap_r = rho_p
    n = np.linspace(0,ncut,ncut+1,dtype=int)
    bes = scs.spherical_jn(n,beta_b)
    bes_1 = scs.spherical_jn(n,beta_b,derivative=True)
    bes_b = scs.spherical_jn(n,beta_pb)
    bes_1_b = scs.spherical_jn(n,beta_pb,derivative=True)       
    neu = scs.spherical_yn(n,beta_b)
    neu_1 = scs.spherical_yn(n,beta_b,derivative=True)
    han = bes + zi * neu
    han_1 = bes_1 + zi * neu_1
    neu_b = scs.spherical_yn(n,beta_pb)
    neu_1_b = scs.spherical_yn(n,beta_pb,derivative=True)
    han_b = bes_b + zi * neu_b
    han_1_b = bes_1_b + zi * neu_1_b        
    bes_p = scs.spherical_jn(n,beta_pb) + th[n] * han_b            
    bes_p_1 = scs.spherical_jn(n,beta_pb,derivative=True) + th[n] * han_1_b
    tz = kap_k / kap_r * bes * bes_p_1 - bes_1 * bes_p
    tn = kap_r / kap_k * han_1 * bes_p - bes_p_1 * han
    tm_2l_hp = kap_r / kap_k * tz / tn
    return tm_2l_hp
   
























