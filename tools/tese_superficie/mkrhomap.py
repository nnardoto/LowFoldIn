#!/usr/bin/env python3
"""mkrhomap.py pairs prefix y0 ny nper_y [x0 z0]: nml do rhomap mode='cube' (x: 48 pts em a, y: a2y/nper_y, z: 60 pts em c)."""
import sys
pairs, pre, y0, ny, nper = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
x0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
z0 = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
A, A2Y, C = 4.7672460267921517, 24.771328996640282 / 6, 5.9335974928231314
NX, NZ = 48, 60
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
# raízes configuráveis (padrão: layout original da campanha)
LFI_TESE = _os.environ.get("LFI_TESE", "/home/nathanael/calc/tese-Te/green_surface/lfi_tese/")
PAO = _os.environ.get("LFI_PAO", _os.path.join(_HERE, "..", "..", "test", "data", "pao", "Te7.0.pao"))
print(f"""&RHOMAP source = 'pairs', pairs_file = '{pairs}', pairs_window = 'occ', spin_factor = 2.0, output_prefix = '{pre}' /
&SPECIES nspec = 1, label = 'Te', pao_file = '{PAO}', basis = 's3p3d2f1', rcut_bohr = 7.0 /
&GRID mode = 'cube', units = 'ang', origin = {x0:.12f}, {y0:.12f}, {z0:.12f},
  du = {A/NX:.15f}, 0.0, 0.0, dv = 0.0, {A2Y/nper:.15f}, 0.0, dw = 0.0, 0.0, {C/NZ:.15f},
  nu = {NX}, nv = {ny}, nw = {NZ} /""")
