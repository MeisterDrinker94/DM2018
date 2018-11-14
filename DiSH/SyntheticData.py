import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def createSynthetic(clusterPoints=300,noisePoints=150,x=0.2,y=0.8):
	"""	
		clusterPoints = number of Points in Clusters
		noisePoints = number of noise points
		x = intersection of y-parallel cluster
		y = intersextion of x-parallel cluster
		creates a synthetic DataSet in 3 dimensional space, with 3 planes
		of points intersecting each other in the same line
	"""
	
	#generate random values
	helper = np.random.rand(clusterPoints,5)
	
	#initialize empty list
	data = []
	
	for ar in helper:
		data.append([x,ar[0],ar[1]])
		data.append([ar[2],y,ar[3]])
		data.append([x,y,ar[4]])
	
	noise = np.random.rand(noisePoints,3)
	
	for ar in noise:
		data.append(ar)
		
	return np.array(data)
		
def main():
	d = createSynthetic(300,150,0.5,0.2)
	
	fig = plt.figure()
	
	ax = fig.add_subplot(111, projection='3d')
	
	ax.scatter3D(d[:,0],d[:,1],d[:,2])
	plt.show()
		
		
if __name__ == '__main__':
	main()	
