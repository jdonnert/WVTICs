# Toycluster

Idealized Cluster Mergers for Gadget (Donnert 2014)

	* Physical model from Donnert 2014
	* Simple OpenMP Tree for neighbour finding with WC6
	* DM distribution function using numerical solution of Eddingtons eq.
	* SPH density sampling & relaxation with Voronoi Tesselations (Diehl+)
	* Magnetic fields from vector potential, Model from Bonafede+ 2010
	* Merger using Comet like or parabular shapes
	* Velocities parametrized using zero-energy orbit
	* optional: substructure population from Giocoli 2010
	* optional: placement of a third halo anywhere in the box
	
Known Issues: 
	* Magnetic field normalisation resolution dependent, due to low accuracy
		SPH rotation operator. DivB=0 constraint pretty bad.
