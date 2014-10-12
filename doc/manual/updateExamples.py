import os

dirGallery = '/home/ratnani/Projects/femiga/branches/0.9.2/gallery/'

files = ['demo/basic_pde_1d', 'demo/basic_pde_2d', 'demo/dirichlet_1d',
         'demo/dirichlet_2d', 'demo/neumann_2d',
         'demo/dirichlet_neumann_2d', 'anisotropicDiffusion',
         'bilaplacian', 'nonlin_ex1_picard', 'nonlin_ex1_newton', 'monge_ampere_picard']

for f in files:
    filename = os.path.join(dirGallery, f)
    os.system('cp '+ filename + '.py'+ ' examples')
    try:
        os.system('cp '+ filename + '.png'+ ' examples')
    except:
        pass


files = ['demo/test_solver_scipy.py', 'demo/test_solver_petsc4py.py', 'demo/test_solver_pyamg.py']
for f in files:
    filename = os.path.join(dirGallery, f)
    os.system('cp '+ filename + ' include/demo')
