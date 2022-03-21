# compiler settings
fortran_compiler = gfortran

# tools and options
delete = rm -rf
makedir = mkdir -p
print = echo

# directories
source_directory = src
module_directory = mod
binary_directory = bin
object_directory = obj
output_directory = out

# extensions
source_extension = .f90
module_extension = .mod
binary_extension = .exe
output_extension = .plt
object_extension = .o

# source, expected object, module and output files
source_files = $(wildcard $(source_directory)/*$(source_extension))
object_files = $(patsubst $(source_directory)/%$(source_extension), \
	$(object_directory)/%$(object_extension), $(source_files))
module_files = $(patsubst $(object_directory)/%$(object_extension), \
	$(module_directory)/%$(module_extension), $(object_files))
output_files = $(wildcard $(output_directory)/*$(output_extension))

# program source and expected object files
source_file = $(wildcard *$(source_extension))
object_file = $(patsubst %$(source_extension), \
	$(object_directory)/%$(object_extension), $(source_file))

# executable
executable = $(patsubst %$(source_extension), \
	$(binary_directory)/%$(binary_extension), $(source_file))

# application execution
all: dirs $(executable)
	$(executable)

# application building
$(executable): $(object_file) $(object_files)
	$(fortran_compiler) -o $@ $^

$(object_file): $(source_file) $(object_files)
	$(fortran_compiler) -J $(module_directory) -o $@ -c $<

$(object_directory)/%$(object_extension): $(source_directory)/%$(source_extension)
	$(fortran_compiler) -J $(module_directory) -o $@ -c $<

# create supple, binary, module and object directories
dirs:
	@$(makedir) $(binary_directory) $(module_directory) $(object_directory) \
		$(output_directory)

# delete generated object and module files
clean:
	@$(delete) $(object_files) $(object_file) $(module_files) $(output_files) \
		$(executable)

# print source, generated object and module files
debug:
	@$(print) source_files: $(source_files)
	@$(print) object_files: $(object_files)
	@$(print) module_files: $(module_files)
	@$(print) output_files: $(output_files)

# declare targets that are not files
.PHONY: all clean debug dirs

# dependent modules
obj/functions.o: obj/constants.o obj/matrix_operations.o
