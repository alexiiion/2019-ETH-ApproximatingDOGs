import svgpathtools.path
import sys,os
sys.path.insert(0, os.getcwd() + "/../../libigl/python")
sys.path.insert(0, os.getcwd() + "/../../libigl/external/nanogui/build/python")
import pyigl as igl
from iglhelpers import *
import svgwrite

def submesh_to_svg(in_path, out_path):
	# Load a mesh in OBJ format
	V = igl.eigen.MatrixXd()
	F = igl.eigen.MatrixXi()
	igl.readOBJ(in_path, V, F)
	# Find the open boundary
	bnd = igl.eigen.MatrixXi()
	igl.boundary_loop(F, bnd)

	V = e2p(V)
	V = 10*V
	bnd_vn = bnd.rows()
	bnd = e2p(bnd)
	#polygon_path = svgpathtools.path.Path(*[svgpathtools.path.Line(V[i][0:2], V[(i + 1) % len(V)][0:2]) for i in range(len(V))])
	#dwg = svgwrite.Drawing(out_path, profile='tiny')
	#dwg = svgwrite.Drawing(out_path, height=1000, width=1000)
	dwg = svgwrite.Drawing(out_path, size=('13cm', '14cm'), viewBox=('-200 -200 200 200')) 
	#Drawing(height='10cm', width='20cm')
	for i in range(0,bnd_vn):
		cur_idx, next_idx = bnd[i], bnd[(i+1)%bnd_vn]
		start_x, start_y = V[cur_idx][0][0:2]
		end_x, end_y = V[next_idx][0][0:2]
		dwg.add(dwg.line((start_x,start_y), (end_x, end_y), stroke=svgwrite.rgb(10, 10, 16, '%')))
	dwg.save()


if __name__ == "__main__":
	if len(sys.argv) == 3:
		folder_path = sys.argv[1]
		model_name = sys.argv[2]
		for filename in os.listdir(sys.argv[1]):
			if (filename.find("__parts") != -1) and filename.endswith("-uv.obj"):
				submesh_n_idx1,submesh_n_idx2 = filename.find("parts")+len("parts"),filename.find("-uv.obj")
				submesh_n = filename[submesh_n_idx1:submesh_n_idx2]
				input_mesh, output_mesh = folder_path+ "//" + filename, folder_path+ "//" + "parts" + submesh_n + ".svg"
				submesh_to_svg(input_mesh, output_mesh)
	else:
		print "Usage: submesh_bnd_to_svg folder_path mesh_name"