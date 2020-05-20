from sfepy.discrete import Problem
from sfepy.base.base import IndexedStruct
from sfepy.postprocess.viewer import Viewer
from utilities import *
import os

acoustics_file = 'data/acoustics.py'
output_solution = 'data/acoustics.vtk'

def resolve_acoustics_problem():
    # Import problem description file
    print_header('Importing problem file')
    pb = Problem.from_conf_file(acoustics_file)

    # Solve the problem
    print_header('Resolution')
    status = IndexedStruct()
    state = pb.solve(status=status)

    # Save its solution
    pb.save_state(output_solution, state)
    os.system('del room.vtk')

def view_acoustics_solution():
    # View the results of the problem resolution
    print_header('View')
    view = Viewer(output_solution)
    view(vector_mode='warp_norm', rel_scaling=2, is_scalar_bar=True, is_wireframe=True)