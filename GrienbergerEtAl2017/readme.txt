I) To install on a local machine:

1) install anaconda python, includes ipython, numpy, matplotlib
https://docs.continuum.io/anaconda/install

add additional ipyparallel package:
conda install ipyparallel

2) download neuron from:
http://www.neuron.yale.edu/ftp/neuron/versions/alpha/nrn-7.4.rel-1390.tar.gz

or install from mercurial
cd ~/neuron
hg clone http://www.neuron.yale.edu/hg/neuron/nrn -r 'Release 7.4'
hg clone http://www.neuron.yale.edu/hg/neuron/iv

following tips from:

http://www.neuron.yale.edu/neuron/download/compilestd_osx

http://www.neuron.yale.edu/phpBB/viewtopic.php?f=4&t=3051#p12584

and making sure to execute these commands before running 'make':

export CFLAGS='-Qunused-arguments'
export CXXFLAGS='-Qunused-arguments'

cd nrn/src/nrnmpi
sh mkdynam.sh

3) Install btmorph from

http://btmorph.readthedocs.org/en/latest/readme.html#installation

4) Make sure ~/neuron/nrnenv includes:

export IDIR=/Applications/NEURON-7.4
export IV=$IDIR/iv
export N=$IDIR/nrn
export CPU=x86_64
export PATH=$IV/$CPU/bin:$N/$CPU/bin:$PATH


and ~/.bash_profile includes:

export HOME=/Users/milsteina
export PATH="$HOME/anaconda/bin:$PATH"
export PATH=$HOME/local/bin:$PATH
source $HOME/neuron/nrnenv
export PATH="/Applications/NEURON-7.4:$PATH"
export LD_LIBRARY_PATH=$HOME/local/lib/
export LD_LIBRARY_PATH="$HOME/anaconda/lib:$LD_LIBRARY_PATH"
export PYTHONHOME="$HOME/anaconda"
export PYTHONPATH="$PYTHONPATH:/Applications/NEURON-7.4/nrn/lib/python"
export PYTHONPATH="$PYTHONPATH:$HOME/python_modules/btmorph"
export PATH="$PYTHONHOME:$PYTHONPATH:$LD_LIBRARY_PATH:$PATH"

5) If error related to libreadline.6.2.dylib :

cd /Applications/NEURON-7.4/nrn/lib/python/neuron/
install_name_tool -change libreadline.6.2.dylib $HOME/anaconda/lib/libreadline.6.2.dylib hoc.so

II) To install on a linux cluster:

1) Download the latest tar.gz from http://www.neuron.yale.edu/ftp/neuron/versions/
e.g.
wget http://www.neuron.yale.edu/ftp/neuron/versions/alpha/nrn-7.4.rel-1324.tar.gz

2) Unpack the file:
tar xvzf nrn-7.4.rel-1324.tar.gz

3) create these subdirectories in your home directory on the cluster: /neuron/nrn

4) install following these directions:
http://www.neuron.yale.edu/neuron/download/compile_linux

from the directory containing the unzipped install files:
./configure --prefix=$HOME/neuron/nrn --without-x --with-paranrn=dynamic --with-nrnpython=dynamic
make
make install

- produces some long pauses with error messages related to a deprecated NumPy API, I think this is a strange interaction with anaconda python2.7, but it doesn't appear to cause any problems.

5) Make sure ~/neuron/nrnenv includes:

export IDIR=$HOME/neuron
export N=$IDIR/nrn
export CPU=x86_64
export PATH=$N/$CPU/bin:$PATH

6) Make sure ~/.bash_profile includes:

export PATH=/usr/local/anaconda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/anaconda/lib
source $HOME/neuron/nrnenv
export PATH=$HOME/neuron/nrn:$PATH
export PYTHONHOME=/usr/local/anaconda
export PYTHONPATH=$PYTHONPATH:$HOME/neuron/nrn/lib/python
export PYTHONPATH=$PYTHONPATH:$HOME/python_modules/btmorph
export PYTHONPATH=$PYTHONPATH:$HOME/CA1Sim
export PATH=$PYTHONHOME:$PYTHONPATH:$LD_LIBRARY_PATH:$PATH

7) mpi4py is required to use the IPython.parallel framework across cores on more than one node, and it is already installed alongside anaconda on the Janelia cluster.
In order to use the Intel MPI implementation on the Janelia cluster, make sure the following lines are included in ~/.bashrc

if [ -f /usr/local/INTEL2016.sh ]; then
  . /usr/local/INTEL2016.sh
fi

export I_MPI_DEBUG=0
export I_MPI_RENDEZVOUS_RDMA_WRITE=1
export I_MPI_DEVICE=rdssm:OpenIB-iwarp
export I_MPI_FALLBACK_DEVICE=0
export I_MPI_USE_DYNAMIC_CONNECTIONS=0

8) In order for the IPython use more than one node on the cluster

A custom ipython profile must be used for the IPython.parallel framework to work with an MPI backend on the Janelia cluster.

ipython profile create --parallel --profile=mpi
edit the file IPYTHONDIR/profile_mpi/ipcluster_config.py
c.IPClusterEngines.engine_launcher_class = 'MPIEngineSetLauncher'
c.IPClusterEngines.work_dir = u'/groups/magee/home/milsteina/CA1Sim/'

edit the file IPYTHONDIR/profile_mpi/ipengine_config.py
c.MPI.use = 'mpi4py'





6) To run my simulations, make sure to execute nrnivmodl in the directory that contains the .mod and .py files.
Make sure to put .swc files in the /morphologies directory and expect to find .pkl and .hdf5 data output files in the
/data directory

Then start an iPython session with: ipython

and execute pieces of the similation with:
run name_of_py_file

build your own cell with

from specify_cells import *
from function_lib import *
from plot_results import *

cell = HocCell()
cell.make_section('soma')

You can also import a morphology from an .swc file:

cell = HocCell('str_with_name_of_swc_file.swc')

BtMorph requires that the compartments with indices 1, 2, and 3 all be of type 1 (soma), according to the standard used by NeuroMorpho.org . My model expects compartment types to be labeled as follows (1 = soma, 2 = axon, 3 = basal, 4 = apical, 5 = trunk, 6 = tuft). Soma and axon compartments are discarded and replaced with a simplified representation of soma and axon.

Ion channel mechanisms and cable parameters for various types of compartments can be inserted with commands like:

cell.modify_mech_param('soma', 'cable', 'Ra', 150.)
cell.modify_mech_param('soma', 'pas', 'g', 0.0002)
cell.modify_mech_param('basal', 'h', 'ghbar', origin='soma')
cell.modify_mech_param(sec_type, 'kap', 'gkabar', origin='soma', slope=3.84e-4, max=0.24)

Once you have a cell appropriately specified, you can save a file containing the mechanisms with:

cell.export_mech_dict('str_with_name_of_pkl_file.pkl'), or with no arg to automatically generate a unique filename with
the data and timestamp.

Then a cell can be instantiated with a mech_dict file:

cell = HocCell('str_with_name_of_swc_file.swc', 'str_with_name_of_pkl_file.pkl')

A quick simulation using adaptive timestep integration can be run with:

sim = QuickSim(duration)
sim.append_rec(cell, cell.tree.root, loc=0.5, description='soma')
sim.run()
sim.plot()

or

sim.export_to_file(name_of_hdf5_object, index_of_similation_in_file)

then calling a plot function on the output file:

plot_superimpose_conditions('str_with_name_of_output_file') # no suffix necessary


7) PyCharm limits the size of the console output buffer. Change the value of idea.cycle.buffer.size in the idea.properties file in the /bin directory of the install package. To change the size of the terminal output buffer, change the registry key terminal.buffer.max.lines.count. Navigate to Help| Find action| Type "Registry"| Find terminal.buffer.max.lines.count.

