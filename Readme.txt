  _____  ____   _____   _____  _____   ____       _ ______ _____ _______ 
 |  __ \|  _ \ / ____| |  __ \|  __ \ / __ \     | |  ____/ ____|__   __|
 | |__) | |_) | (___   | |__) | |__) | |  | |    | | |__ | |       | |   
 |  ___/|  _ < \___ \  |  ___/|  _  /| |  | |_   | |  __|| |       | |   
 | |    | |_) |____) | | |    | | \ \| |__| | |__| | |___| |____   | |   
 |_|    |____/|_____/  |_|    |_|  \_\\____/ \____/|______\_____|  |_|   
==========================================================================

Project "rules"
------------------

-> The master branch ALWAYS works
	- compiles
	- gives valid results

-> Create a new branch on the online repo for each feature you are adding
	- Work on that branch until it's working
	- When finished, merge the master into your own branch and make sure everything works correctly.
	- When you have merged the master into your branch and it works correctly, merge it back into the master

-> To tell others when you are working on something, to avoid work duplication

-> Try to keep nice code style/informative comments

-> Push results often. Try to make the last commit of every push have information on the changes that have been made 
	(Ideally all the commit messages should have useful information, but that's hard to keep)

Structure of code (Mainly fluid simulation for now):
------------------
The fluid simulation code should be encapsulated in it's own function. The function takes the data structures corresponding
to step n, and returns the data structures corresponding to step n+1.

The SPH algorithm can be structures as a number of steps, applied to all particles. If each step is encapsulated into a function,
it will be easier to implement different versions of these.

The function may look something like this

void compute_sph_step( data *in, data *out )
//data in/out may actually mean multiple data structures in each place
//The pointers for in and out may be the same, if we modify directly the structures without copying
{
	step_1(...);
	step_2(...);
	..
	step_m(...);
}

This structure can be very useful for working in multiple parts of the implementation without affecting other's code.
Also good, because this steps can be implemented in many different ways. Should allow to interchange them easily if we
implement more than one.

With this, we can implement a render step or store the info at each step in a simple way:

//Pseudocode of structure, memory management stuff is ignored
void render()
{
	data *in,*out;
	set_starting_conditions(in);
	for( i:(1..nsteps) ) {
		compute_sph_step(in,out);
		render_particles(out); //For debug rendering in opengl
		//store_particle_data(out); //For rendering in external source
		in = out;
	}
}


