import numpy as np
import math 

filename = "single.xyz"
xyz = open(filename, "r")	
back = open('single_out.xyz', 'w')

Geometry_Data = xyz.readline() + xyz.readline()
back.write(Geometry_Data)

OP_f = []
elements = []
 
a = 0

c = h = n = zn = 0

mass_carbon = 12.011
mass_zinc = 65.38
mass_hydrogen = 1.008
mass_nitrogen = 14.007

C_x = H_x = N_x = Zn_x = C_y = H_y = N_y = Zn_y = C_z = H_z = N_z = Zn_z = 0

Data = xyz.readlines()		

for line in Data:
	mid = line
	line = mid.split()
	numbers = line[1:4]
 
 	convert_float = [float(i) for i in numbers]
	numbers = np.array(numbers)		

	vector_x = float(numbers[0]) 
	vector_y = float(numbers[1]) 
	vector_z = float(numbers[2])
		
	if line[0] == "C":
		c+=1
		CM_Carbon_x = vector_x*mass_carbon # calculates center of mass
		CM_Carbon_y = vector_y*mass_carbon
		CM_Carbon_z = vector_z*mass_carbon
		C_x += CM_Carbon_x # sums center of mass for each atom
		C_y += CM_Carbon_y
		C_z += CM_Carbon_z
						
	if line[0] == "H":
		h+=1
		CM_Hydrogen_x = vector_x*mass_hydrogen
		CM_Hydrogen_y = vector_y*mass_hydrogen
		CM_Hydrogen_z = vector_z*mass_hydrogen
                H_x += CM_Hydrogen_x
		H_y += CM_Hydrogen_y
		H_z += CM_Hydrogen_z
				

	if line[0] == "N":
		n+=1
		CM_Nitrogen_x = vector_x*mass_nitrogen
                CM_Nitrogen_y = vector_y*mass_nitrogen
		CM_Nitrogen_z = vector_z*mass_nitrogen
		N_x += CM_Nitrogen_x
		N_y += CM_Nitrogen_y
		N_z += CM_Nitrogen_z
		
	if line[0] == "Zn":
		zn +=1
		CM_Zinc_x = vector_x*mass_zinc
		CM_Zinc_y = vector_y*mass_zinc
		CM_Zinc_z = vector_z*mass_zinc
                Zn_x += CM_Zinc_x
		Zn_y += CM_Zinc_y
		Zn_z += CM_Zinc_z
		
Total_Mass = mass_carbon*c + mass_zinc*zn + mass_hydrogen*h + mass_nitrogen*n
Center_Mass_x = (C_x + H_x + N_x + Zn_x)/Total_Mass # Radius of center of mass calculation for each coordinate
Center_Mass_y = (C_y + H_y + N_y + Zn_y)/Total_Mass		
Center_Mass_z = (C_z + H_z + N_z + Zn_z)/Total_Mass	
print "Center Mass x:", Center_Mass_x, "Center Mass y:", Center_Mass_y, "Center Mass z: ", Center_Mass_z
print c, n, zn, h

while a < 1.8:
	for l in Data:
	
		st = l
		l = st.split()

		x = math.pi*a
		s = math.sin(x)
		c = math.cos(x)

		elements.append(l[0])

		fix = l[1:4]	#take only numerical entries

		final = [float(i) for i in fix] #convert string to float
		final = np.array(final)	# make an array of the floating point number

		y = 0.152593 	# 1.9 - 1.462967071 (angle that seems to work) 
				# angle of rotation for the molecule at the origin
                si = math.sin(y)
                co = math.cos(y)
		
		final[0] = final[0] - Center_Mass_x # move molecule to origin
		final[1] = final[1] - Center_Mass_y
		final[2] = final[2] - Center_Mass_z
		B1=np.matrix([[co, -si, 0], [si, co, 0], [0, 0, 1]]) #rotates vector CCW
		OP_R = np.dot(B1, final)
		
		OP_R = np.array(OP_R)
		OP_R[0][0] = OP_R[0][0] + float(21.08239343)  # add coordinates back to molecule for circular shape
                OP_R[0][1] = OP_R[0][1] + float(4.627052911) 
                #OP_R[0][2] = OP_R[0][2]  
		OP_1D = OP_R.flatten()	
		
		final[0] = final[0] + Center_Mass_x # return values to array for next iteration
                final[1] = final[1] + Center_Mass_y 
                final[2] = final[2] + Center_Mass_z 
		 
		B2=np.matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]]) # transformation matrix
		OP=np.dot(B2, OP_1D)	# perform dot product on the matrices 	
		
 		OP_f.append(OP)         #creates an entire matrix to include all lines 
	
		
	a += float(0.2)

OP_f = np.array(OP_f)		# converts to numpy array

OP_i = OP_f.flatten()		# creates a 1-D array  

OP_f = [OP_i[i:i+3] for i in xrange(0,len(OP_i),3)] 	# 'shapes' the array such that there are 3 columns
 	
Elements = np.array(elements)

OP_f_text = np.savetxt('out.tmp', OP_f) 	# reads array and saves into a file as a string
OP_f_infile = open('out.tmp','r')	# open file to write string array onto
OP_f_text = OP_f_infile.read()	

line_by_line = OP_f_text.split('\n')		# rearrange text file such that coordinates are readable line by line
line_by_line = np.array(line_by_line)		# convert to array

zipped = zip(Elements, line_by_line)		# combine the element array and coordinate array in pairs 			
zipped2 = [i + ' ' + j for i,j in zipped]	# rearrange so that each pair has a line
zipped2 = np.array(zipped2)			# convert to array

final = np.savetxt('output.tmp', zipped2, fmt = '%s')	# convert to a text file so that it is readable
infile = open('output.tmp','r')
final = infile.read()

back.write(final)		# write content onto output file

back.close()
xyz.close()




