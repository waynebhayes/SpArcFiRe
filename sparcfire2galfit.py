import sys
import math
from spiralgalaxy import SpiralGalaxy

class Sparcfire2Galfit:

	def __init__(self,sg):
		self.spiral_galaxy = sg

	def output_galfit(self):
		original_stdout = sys.stdout
		f = open(self.spiral_galaxy.name + "_S2G.txt",'w')
		sys.stdout = f

		self.control_parameters()
		self.bulge()
		# self.arcs()
		self.sky()

		sys.stdout = original_stdout
		f.close()

	def control_parameters(self):
		print("""
A) gal.fits            # Input data image (FITS file)
B) imgblock.fits       # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D) psf.fits   #        # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1    93   1    93   # Image region to fit (xmin xmax ymin ymax)
I) 100    100          # Size of the convolution box (x y)
J) 26.563              # Magnitude photometric zeropoint 
K) 0.038  0.038        # Plate scale (dx dy)    [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps
		""")


	def bulge(self):
		print("#Object Number 1")
		functions = ["object_type","position","integrated_magnitude","half_light_radius",
					"sersic_index","axis_ratio","position_angle","rotation_function","R1",
					"R2","R3","R4","R9","R10","B3","Z"]
		for f in functions:
			print(eval("self." + f + "()"))

	def arcs(self):
		for arc in spiral_galaxy.arcs:
			pass

	def sky(self):
		print("""
# Object number: 2
0) sky                    #  object type
1) 5.49		1          #  sky background at center of fitting region [ADUs]
2) 0.0000      0          #  dsky/dx (sky gradient in x)
3) 0.0000      0          #  dsky/dy (sky gradient in y)
Z) 0                      #  output option (0 = resid., 1 = Don't subtract)
		""")

	def object_type(self):
		return "0)\tsersic"

	def position(self):
		return "1)\t48.5180\t51.2800\t1\t1"

	def integrated_magnitude(self):
		# m = u - 5log(a) - 2.5log[(2n)!(pi)]
		u = 0 # surface brightness
		n = 2 # sersic index
		a = 0 # scale length
		m = u - (5 * math.log(a,10)) - (2.5 * log(math.factorial(2 * n) * math.pi))
		return "3)\tstr(m)\t1"

	def half_light_radius(self):
		return "4)\t4.5160\t1"

	def sersic_index(self):
		return "5)\t8.2490\t1"

	def axis_ratio(self):
		return "9)\t0.7570\t1"

	def position_angle(self):
		return "10)\t60.3690\t1"

	def rotation_function(self):
		return "R0)\tpowerlaw"

	def R1(self):
		return "R1)\t30\t1"

	def R2(self):
		return "R2)\t100\t1"

	def R3(self):
		return "R3)\t275\t1"

	def R4(self):
		return "R4)\t0.5\t1"

	def R9(self):
		return "R9)\t0.5\t1"

	def R10(self):
		return "R10)\t0.5\t1"

	def B3(self):
		return "B3)\t0.03\t1"

	def Z(self):
		return "Z)\t0"
