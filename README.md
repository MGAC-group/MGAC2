# The Modified Genetic Algorithm for Crystals (MGAC)

This code is a modernization of the MGAC code written by the Facelli group. The
code has been contributed by a number of authors over the years:

- Victor E. Bazterra
- Gabriel I. Pagola
- Albert M.Lund
- Anita M. Orendt 
- Marta B. Ferraro
- Julio C. Facelli

Users of this code should cite the following publications:

```
Victor E. Bazterra, Marta B. Ferraro, and Julio C. Facelli. The Journal of Chemical Physics 116, 5984 (2002); doi: http://dx.doi.org/10.1063/1.1458547

Lund AM, Pagola GI, Orendt AM, Ferraro MB, Facelli JC. Crystal Structure Prediction from First Principles: The Crystal Structures of Glycine. Chemical physics letters. 2015;626:20-24. doi:10.1016/j.cplett.2015.03.015.
```

Scientific or collaborative inquiries should be directed to Julio Facellio (julio.facelli@utah.edu)

For assistance with the code or other technical aspects please contact Albert Lund (a.m.lund@utah.edu)



================================================================================
## INSTALL INSTRUCTIONS:

Requires:

- cmake 2.8.12 or greater
- mpich2 3.1.4 or greater
- gcc 4.9.2 or greater (C++11 support)
- QE 5.0.2 (other versions NOT tested)
- objcopy (should be standard on linux)

Intel compilers/MPI optional but not required
- Intel compilers v2015.1.133 or greater
- IMPI 5.0.1.035 or greater

On generic systems without Intel compilers:
```
mkdir mgac2-build
cd mgac2-build
cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc ../mgac2
make
```
MAKE SURE TO SOURCE MPICH2 VARS!

On CHPC systems with Intel compilers:

```
module load cmake intel impi gcc/4.9.2
mkdir mgac2-build
cd mgac2-build
cmake -D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc ../mgac2
make
```

MODULE LOAD ORDER IS IMPORTANT! gcc/4.9.2 must be loaded last

To run make sure that the gcc/4.9.2 libraries are accessible. If using mpi 3.1.4 or greater mgac can use either mpich2 or impi to run. 

================================================================================
## RUNNING NOTES:

To see the help message:
`mgac.x -h`
`mgac.x --help`

A simple run is executed by:
`mgac.x -i inputfile`

A restart is the same, with a reference to the sql file:
`mgac.x -i inputfile -r restart.sq3`

A list of spacegroups can be listed by:
`mgac.x -l `

To convert an output sqlite file to CIF format. This automatically sorts structures:
`mgac.x -i output.sq3 -c file.cif`

The number of best structures for CIF conversion defaults to 100. Use the -s flag to specify a different number:
`mgac.x -i output.sq3 -c file.cif -s 200`

A different table can be selected from the default "structs" table using -b:
`mgac.x -i output.sq3 -c file.cif -s 100 -b precluster`

An input template can be quickly generated from a CIF using -t:
`mgac.x -i input.cif -t template.xml`

The input template still requires all other parameters to be manually added! 

The plane for molecules must be specified using -p in comma delimited form using the atom labels specified in the CIF:
`mgac.x -i input.cif -t template.xml -p C1,C2,N1`

NOTE: the template generation does not include bond/angle/limitations or dihedrals. Those must be specified manually still.

================================================================================
## INPUT FORMAT

The MGAC2 format uses a simplified XML format over the MGAC1 format. An example is given below:

```
<mgac>
        <crystal interdist="0.0" contactdist="0.3" intradist="1.4" maxvol="100.0" minvol="-30.0" name="HISTAN">
                <cell>
                        <stoichiometry mol="gly-1" count="1"/>
                </cell>
                <molecule name="gly-1" plane="C1 C2 N1" expectvol="156.6">
                        <atom elem="C" title="C1" x="2.35600" y="-0.72600" z="-0.12700"/>
                        <atom elem="C" title="C2" x="1.35200" y="1.13200" z="-0.05100"/>
                        <atom elem="C" title="C3" x="0.41200" y="0.21500" z="0.27100"/>
                        <atom elem="C" title="C4" x="-1.04300" y="0.33800" z="0.60200"/>
                        <atom elem="C" title="C5" x="-1.97500" y="-0.21000" z="-0.49200"/>
                        <atom elem="H" title="H1" x="0.69400" y="-1.88700" z="0.40700"/>
                        <atom elem="H" title="H2" x="3.08900" y="-1.49900" z="-0.23500"/>
                        <atom elem="H" title="H3" x="1.23300" y="2.19400" z="-0.11900"/>
                        <atom elem="H" title="H4" x="-1.26300" y="1.38800" z="0.76500"/>
                        <atom elem="H" title="H5" x="-1.25700" y="-0.17100" z="1.54200"/>
                        <atom elem="H" title="H6" x="-1.77200" y="0.30800" z="-1.42300"/>
                        <atom elem="H" title="H7" x="-1.75600" y="-1.25800" z="-0.67200"/>
                        <atom elem="H" title="H8" x="-3.63300" y="-0.57000" z="0.64100"/>
                        <atom elem="H" title="H9" x="-3.66300" y="0.86300" z="-0.09900"/>
                        <atom elem="N" title="N1" x="1.07600" y="-0.99000" z="0.21700"/>
                        <atom elem="N" title="N2" x="2.56200" y="0.53300" z="-0.29800"/>
                        <atom elem="N" title="N3" x="-3.39300" y="-0.09600" z="-0.20800"/>

                        <dihedral title="C3-C4-C5-N3" angle="C3 C4 C5 N3" update="N3 H8 H9 H6 H7" />
                        <dihedral title="C2-C3-C4-C5" angle="C2 C3 C4 C5" update="C5 N3 H8 H9 H6 H7 H4 H5" />
                        <dihedral title="C4-C5-N3-H8" angle="C4 C5 N3 H8" update="H8 H9" />

                </molecule>
        </crystal>
        <qe>
                <param calculation="vc-relax" restart_mode="from_scratch" 
                                tstress=".true." tprnfor=".true." 
                                nstep="150"
                                qepath="pw.x"
                                prefix="temp"
                                mpirunpath="mpirun -bootstrap ssh -rmk user "
                                preamble="module load intel/2015.1.133 impi/5.0.1.035 qe/5.0.2"
                                pseudo_dir="/uufs/chpc.utah.edu/sys/installdir/qe/espresso-5.0.2/pseudo/" 
                                outdir="/scratch/general/lustre/u0403692/mgac" 
                                wf_collect="true" verbosity="low" 
                                etot_conv_thr="1.0D-3" forc_conv_thr="1.0D-2" press_conv_thr="0.5D0" 
                                ecutwfc="55" ecutrho="550" spline_ps=".true." 
                                conv_thr="1.D-7" cell_dynamics="bfgs" 
                                k_points="automatic" k_point_spec="2 2 2 1 1 1" 
                                restart_limit="3" scftimeout="7200" />


                <pseudo elem="H" mass="1.000" name="H.pbe-rrkjus.UPF" />
                <pseudo elem="N" mass="15.000" name="N.pbe-rrkjus.UPF" />
                <pseudo elem="C" mass="12.000" name="C.pbe-rrkjus.UPF" />
        </qe>
        <run replacement="2.5" mutation="0.001" generations="100" popsize="30" seed="-1" const_scale="1.0" lin_scale="1.0" exp_scale="0.0" calcmethod="qe" mode="stepwise" type="clustered" selector="roulette" outputmode="sql" outputfile="histan" spacemode="limited" group="4" symmlimit="4" binlimit="10" downlimit="4" precompute="50" clusterdiff="0.3" clustersize="2.0" angstep="2.5"/>
</mgac>
```

All keywords must be encapsulated in an mgac tag. 

Explanation of keywords:

- `crystal::interdist` : When performing fitcell, interdist is the buffer distance between the Van der Waals shells (in Angstroms) between atoms in different molecules

- `crystal::intradist` : When performing fitcell, the buffer distance between Van der Waals shells (in Angstroms) between atoms in the same molecule

- `crystal::maxvol, crystal::minvol` : The maximum and minimum volume percentages for the volume filter.

- `crystal::name` : The name of the crystal.

- `crystal::cell` : Contains informations about the cell. For input this is mainly the stoichiometry, which is currently required but not well implemented for more than one molecule. 

- `crystal::cell::stoichiometry::mol` : The name of the molecule to be used for this stoichimetry. THIS MUST MATCH THE NAME OF AT LEAST ONE MOLECULE!

- `crystal::cell::stoichiometry::count` : The count of the specified molecule.

- `crystal::molecule::name` : The name of the molecule.

- `crystal::molecule::plane` : The space-delimited list of atoms that form a rigid plane in the molecule. The names must match the title of at least three atoms in the atomlist. 

- `crystal::molecule::expectvol` : The estimated volume of the ideal crystal structure. Use arh to obtain this valume (in angstroms^3).

- `crystal::molecule::atom::elem` : The element of the atom. Possible candidates are C,H,N,O,P,Cl,F,S,Br. 

- `crystal::molecule::atom::title` : The label of this atom. ALL ATOM TITLES MUST BE UNIQUE.

- `crystal::molecule::x,y,z` : The x,y, and z coordinates of the atom in cartesian space. Must be given in ANGSTROMS.

- `crystal::molecule::dihedral::title` : The name of a dihedral angle specification. 

- `crystal::molecule::dihedral::angle` : The space-delimited list of atoms that form the dihedral angle. The names of the atoms must form a dihedral angle through bond connectivity and must match atom titles. 

- `crystal::molecule::dihedral::update` : The space-delimited list of atoms that will be rotated for the dihedral. This includes at least one edge atom in the actual dihedral tuple. All atoms must match an atom title. 

- `crystal::molecule::dihedral::min, max` : Minimum and maximum rotation of the dihedral angle, in degrees.

- `qe::param` : A list of QE parameters to be parsed. Many parameters are the same as specified in QE. See the QE website for more info. Only some paramters are discussed here

- `qe::param::nstep` : Number of vc-relax steps before QE exits. 

- `qe::param::qepath` : Full path to QE program. 

- `qe::param::mpirunpath` : Path and arguments to run MPI. 

- `qe::param::preamble` : Command(s) to be run prior to starting QE (e.g. module load or setup script).

- `qe::param::pseudo_dir` : Directory containing psuedopotentials for QE. 

- `qe::param::restart_limit` : The number of times MGAC will restart QE before stopping permanently. A good default is 3.

- `qe::param::scftimeout` : The time limit (in seconds) for a QE run. QE will be forcefully terminated after the time limit. NOTE: If QE is forecfully terminated too frequently, MPI communication may break and nodes will be eliminated from the MGAC run.

- `qe::param::pseudo::elem` : The element of the pseudopotential.

- `qe::param::pseudo::mass` : Mass in AMU of the pseudopotential atom. 

- `qe::param::pseudo::name` : The name of the pseudopotential file.

- `run::replacement` : The fraction replacement for partial cross. E.g., for a population of 50 and replacement of 2.5, 250 structures will be generated at each generation. 

- `run::mutation` : The probability that an individual structure will be mutated. Genes are mutated within that structure at a rate of 0.15. 

- `run::generations` : Number of generations to run over. 

- `run::popsize` : Base population size of run.

- `run::seed` : Seed for the random number generator. If the seed is -1 than a random seed will be chosen. 

- `run::const_scale, lin_scale, exp_scale` : The constant, linear, and exponential scaling factors for the fitness scaling step. The scaling follows the equation scale(score) = c + l*lin_score + e*exp_score where c,l,and e are the scaling factors, respectively, and the score is derived directly from the energy of the structure. The lin_score is is the normalized difference with the minimum score: 1 - [(score - min)/(max/min)]. The exponential_score is derived from the lin_score via E^(-0.5*10.0*lin_score).

- `run::calcmethod` : The calculation method of the system. "qe" spceifies the use of QE, whereas "fitcell" runs using fitcell only.  

- `run::mode` : Can be "stepwise" or "steadystate". Stepwise is the standard generation based GA, steadystate is unimplemented. 

- `run::type` : Can be elitism, fullcross, finaleval, or clustered. Elitism follows the same method as MGAC1. Fullcross uses a fullcrossing method, but can be unstable. Finaleval takes a population as a restart and performs a final QE evaluation on the population. Clustered is the preferred method and eliminates duplicate structures from the structure generation. 

- `run::selector` : The type of structure selector. Only "roulette" is implemented. 

- `run::outputmode` : Only "sql" is implemented. "xml" no longer works. 

- `run::outpufile` : The prefix to use for the output file. Outputfiles will be of the form "outputfile.type", e.g. "histan.sq3".

- `run::spacemode` : Can be single,limited, or full. Single mode selects a single spacegroup and requires the additional use of the "group" keyword. Limited only uses a subset of low complexity spacegroups and is preferred. Full mode uses all spacegroups (NOTE: the last spacegroup is not implemented, so this might cause crashes).

- `run::group` : Specifies the spacegroup number of a single spacegroup run. Use ./mgac.x -l to see a list of spacegroups. 

- `run::symmlimit` : Limits the number of symmetry operations in the available spacegroups. For example, setting a symmlimit of 4 will limit the available spacegroups to spacegroups with 4 or less symmetry operations. 

- `run::binlimit` : In a multiple spacegroup run, spacegroups are sorted into individual spacegroups. When crossing occurs, only the "binlimit" lowest energy spacegroups will be crossed. The energy of a spacegroup is determined from the lowest energy structure in that spacegroup. 

- `run::downlimit` : The number of nodes that can crash before MGAC will terminate. If a node loses communication over MPI it is considered crashed. 

- `run::precompute` : The number of precompute steps that will be performed. A precompute step clusters structures and obtains a maximally diverse set of structures. The number of precompute steps required is dependent on the number of genes, clustersize, and clusterdiff. 

- `run::clusterdiff` : The genetic distance (in fractions of 1.0) between two structures used in clustering. A cluster difference of 0.3 means that 30% of the genes in a structure cannot be shared with any other structure in the clustered populations. 

- `run::clustersize` : When clustering, the gene differences are normalized. Clustersize specifies the effective step size of the normalization. 

- `run::angstep, dihstep, rotstep, ratiostep, posstep` : The step sizes of cell angle, dihedral angle, molecular rotation angles, cell ratios, and fractional positions, respectively, used in the clustering algorithm.
 

================================================================================
## DEVELOPMENT NOTES:

A hook for steady state methods can be found in gasp2.cpp around line 1304, This is where it should likely be implemented, although further changes to the GASPcontrol serverprogram will be needed.



Outlier filter testing: There is a function designed to test if a structure energy is way out of range. (GASP2pop.cpp, line 1496). It isn't working right. 

