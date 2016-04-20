class SpiralGalaxy:

	def __init__(self,name,bar,dar,bab,dmaa,tna):
		self.name = name
		self.bulge_axis_ratio = bar
		self.bulge_avg_brightness = bab
		self.disk_axis_ratio = dar
		self.disk_maj_axis_angle = dmaa
		self.total_num_arcs = tna
		self.arcs = list() #list of spiral arcs

	def create_arc(self,pa,rs,re,al,np,mi):
		self.arcs().append(SpiralArc(pa,rs,re,al,np,mi))

class SpiralArc:

	def __init__(self,pa,rs,re,al,np,mi):
		self.pitch_angle = pa
		self.r_start = rs
		self.r_end = re
		self.arc_length = al
		self.num_pixels = np
		self.mean_intensity = mi