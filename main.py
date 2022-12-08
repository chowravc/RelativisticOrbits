## Importing useful packages
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as co

## Object to create system of relativistic
## body orbiting other mass (M >> m)
class system():

	# Constructor for system
	def __init__(self, e, l, M=u.Msun, dimensions=False):

		# If input values have dimensions
		if dimensions:

			# Set energy and angular momentum
			self.e = (e/(co.c**2)).si
			self.l = l/(co.G*M/co.c)

		# If input values are dimensionless
		else:

			# Set energy and angular momentum
			self.e = e
			self.l = l

		# Precession bias term
		self.eta = np.sqrt(self.l**2 - 1)/self.l

		# Temp variables for eccentricity
		temp1 = (self.e*l)**2 + 1 + 2*self.e*(self.l**2)
		temp2 = (self.e+1)**2

		# Eccentricity
		self.ecc = np.sqrt(temp1/temp2)

		# Delete temp variables
		del temp1; del temp2

		# Radius amplitude
		self.r0 = (self.l**2 - 1)/(self.e+1)

	# Find radii given angle
	def findRs(self, theta):

		# Return radius
		return self.r0/(1 + self.ecc*np.cos(self.eta*theta))

	# Display orbit
	def display(self, thetaI, thetaF):

		# Create list of angles
		thetas = np.linspace(thetaI, thetaF, num=int(1e3*(thetaF-thetaI)/(2*np.pi)))

		# Find list of radii
		rs = self.findRs(thetas)

		# Find x- and y-positions
		xs = rs*np.cos(thetas)
		ys = rs*np.sin(thetas)

		# Plot orbit
		plt.plot(xs, ys)

		plt.gca().set_aspect('equal', adjustable='box')

		# Plotting center
		plt.axhline(0, linestyle='--', color='r')
		plt.axvline(0, linestyle='--', color='r')
		plt.show()

	# Show orbital parameters
	def params(self):

		print(f'e: {self.e:.4g}')
		print(f'l: {self.l:.4g}')
		print(f'ecc: {self.ecc:.4g}')
		print(f'eta: {self.eta:.4g}')
		print(f'r0: {self.r0:.4g}')

## Main functioning of script
def main():

	# Create system
	s1 = system(-1.9, 2.1)

	# Display system parameters
	s1.params()

	# Display orbit for 30 orbits
	s1.display(0, 30*np.pi)

## Called when script is run
if __name__ == '__main__':

	# Run main
	main()