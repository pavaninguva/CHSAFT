from Thermo import PCSAFT
import numpy as np 

r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-3],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
print('da_hc_dx: =======')
print('{:18s}'.format('zeta0_dx: '), r.zetan_dx(0))
print('{:18s}'.format('zeta1_dx: '), r.zetan_dx(1))
print('{:18s}'.format('zeta2_dx: '), r.zetan_dx(2))
print('{:18s}'.format('zeta3_dx: '), r.zetan_dx(3))
print('{:18s}'.format('a_hs: '), r.a_hs())
print('{:18s}'.format('da_hs_dx: '), r.da_hs_dx())
print('{:18s}'.format('dg_hs_dx: '), r.dg_hs_dx())
print('{:18s}'.format('da_HC_dx: '), r.da_HC_dx())
print('da_disp_dx: =======')
print('{:18s}'.format('m2eo3_dx: '), r.m2eno3_dx(1))
print('{:18s}'.format('m2e2o3_dx: '), r.m2eno3_dx(2))
print('{:18s}'.format('I1_dx: '), r.I1_dx())
print('{:18s}'.format('I2_dx: '), r.I2_dx())
print('{:18s}'.format('C1_dx: '), r.C1_dx())
print('{:18s}'.format('da_disp_dx: '), r.da_disp_dx())
print('mu_res: ======')
print('{:18s}'.format('mu_res: '), r.mu_res())
print('{:18s}'.format('a_res: '), r.a_res())
print('{:18s}'.format('Z: '), r.Z())