from sedtools import *
f1 = Flux(5.0,1.0,'test_data/SubB.res')
f2 = Flux(10.0,1.0,'test_data/UKIRTJ.res')
fluxlist = [f1,f2]
redshift = 1.1
name = "test_gal"
preamble =  ['#PBS -l nodes=1','#PBS -l walltime=96:00:00', 
    '#PBS -q lionxf-yuexing','#PBS -j oe', 'cd $PBS_O_WORKDIR',
    'echo " "' , 'echo "Job started on `hostname` at `date`"']
postamble = ['echo "Job Ended at `date`"', 'echo " "']
sedparaminfo = {}
sedparaminfo['Mas']=[16.,9.,35.,0.5]
sedparaminfo['Age']=[16.,12.,23.3,0.5]
sedparaminfo['Tau']=[-0.100,-5.,15.,0.0]
sedparaminfo['EBV']=[0.1,0.0,1.0,0.04]
sedparaminfo['Met']=[-0.7,-2.29,0.45,0.0]
sedparaminfo['Red']=[redshift,0.0,1.5,0.0]

testgal = Galaxy(name, fluxlist, redshift)
galmc = GalMC('test_data/o2default.ini','test_data/empty.depfile')
pbs = PBSBackend(preamble,postamble)

#    def submit(self, backend, galaxy, numchains, sedparaminfo):
galmc.submit(pbs,testgal,4,sedparaminfo)
