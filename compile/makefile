LDLIBS= -I/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/include -L/opt/intel/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64 \
        -mkl
#_lapack -lmkl_em64t -lguide -lpthread -Vaxlib

FC = ifort -c                   #compile
LK = ifort -o agile_idsa        #link
SF = -mcmodel=medium -i-dynamic -no-ipo -r8 -xW -traceback    #compile flags
OF = $(SF) -O2                   #optimization flags
DF = $(SF) -g -O0 -debug semantic_stepping -debug variable_locations -fpstkchk -C -ftrapuv -warn interfaces #debug flags
CF = $(OF)

agile = ../agile/
driver = ../driver/
idsa = ../idsa/

OBJS =\
  adaptive_grid_module.o\
  agile_module.o\
  agile_parameter_module.o\
  agile_start.o\
  atvdflux.o\
  btvdflux.o\
  delept_module.o\
  driver.o\
  egroup_module.o\
  eosinterf.o\
  eos_module.o\
  equations_relax.o\
  equations_step.o\
  etafit.o\
  fermi.o\
  frameint_module.o\
  framerate.o\
  grid_initialisation.o\
  input_file_module.o\
  ionscatt.o\
  meanfreepath.o\
  mu_tau_leakage_module.o\
  newton_module.o\
  nprates.o\
  npscatt.o\
  nuprox.o\
  progenitor.o\
  residuum.o\
  s15s7b2.o\
  sparse_matrix_module.o\
  species_module.o\
  spectrum_module.o\
  state_vector_module.o\
  thermo_module.o\
  units_module.o\
  uydot.o\
  uyupdate.o\
  wrnuprox.o\

agile_idsa:\
  $(OBJS); $(LK) $(CF) $(OBJS) $(LDLIBS)

adaptive_grid_module.o:\
  agile_parameter_module.o\
  input_file_module.o\
  state_vector_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

agile_module.o:\
  agile_parameter_module.o\
  delept_module.o\
  egroup_module.o\
  eos_module.o\
  input_file_module.o\
  newton_module.o\
  species_module.o\
  spectrum_module.o\
  state_vector_module.o\
  units_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)
 
agile_parameter_module.o:\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

agile_start.o:\
  adaptive_grid_module.o\
  input_file_module.o\
  newton_module.o\
  state_vector_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

atvdflux.o:\
  agile_parameter_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

btvdflux.o:\
  agile_parameter_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

delept_module.o:\
  eos_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

driver.o:\
  agile_module.o\
  eos_module.o\
  spectrum_module.o\
  state_vector_module.o\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

egroup_module.o:\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

eosinterf.o:\
  eos_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

eos_module.o:\
  units_module.o\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

equations_relax.o:\
  adaptive_grid_module.o\
  agile_parameter_module.o\
  state_vector_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

equations_step.o:\
  agile_module.o\
  agile_parameter_module.o\
  state_vector_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

etafit.o:\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

fermi.o:\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

frameint_module.o:\
  egroup_module.o\
  species_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

framerate.o:\
  egroup_module.o\
  species_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

grid_initialisation.o:\
  agile_parameter_module.o\
  input_file_module.o\
  state_vector_module.o\
  adaptive_grid_module.o\
  newton_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

input_file_module.o:\
  agile_parameter_module.o\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

ionscatt.o:\
  egroup_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

meanfreepath.o:\
  egroup_module.o\
  eos_module.o\
  species_module.o\
  spectrum_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

mu_tau_leakage_module.o:\
  egroup_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

newton_module.o:\
  adaptive_grid_module.o\
  agile_parameter_module.o\
  input_file_module.o\
  sparse_matrix_module.o\
  state_vector_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

nprates.o:\
  egroup_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

npscatt.o:\
  egroup_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

nuprox.o:\
  egroup_module.o\
  input_file_module.o\
  mu_tau_leakage_module.o\
  species_module.o\
  spectrum_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

progenitor.o:\
  input_file_module.o\
  state_vector_module.o\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

residuum.o:\
  egroup_module.o\
  frameint_module.o\
  species_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

s15s7b2.o:\
  agile_parameter_module.o\
  input_file_module.o\
  units_module.o\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

sparse_matrix_module.o:\
  agile_parameter_module.o\
  input_file_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

species_module.o:\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

spectrum_module.o:\
  egroup_module.o\
  input_file_module.o\
  species_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

state_vector_module.o:\
  agile_parameter_module.o\
  input_file_module.o\
  thermo_module.o\
  units_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

thermo_module.o:\
  agile_parameter_module.o\
  eos_module.o\
  input_file_module.o\
  units_module.o\
  $(agile)$(@:.o=.f); $(FC) $(CF) $(agile)$(@:.o=.f)

units_module.o:\
  $(driver)$(@:.o=.f); $(FC) $(CF) $(driver)$(@:.o=.f)

uydot.o:\
  egroup_module.o\
  species_module.o\
  spectrum_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

uyupdate.o:\
  egroup_module.o\
  frameint_module.o\
  species_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

wrnuprox.o:\
  agile_module.o\
  egroup_module.o\
  input_file_module.o\
  species_module.o\
  spectrum_module.o\
  units_module.o\
  $(idsa)$(@:.o=.f); $(FC) $(CF) $(idsa)$(@:.o=.f)

