

* Copy ssbn.h and ssbn.cpp in directory models/

* Register the new module in models/modelsmodule.cpp
	# include "ssbn.h"	

	kernel().model_manager.register_node_model< ssbn >( "ssbn" );

* Add the model to the list of models which will be compiled for the new version of nest 
	In models/CMakeLists.txt
	set (models_sources
	ac_generator.h ac_generator.cpp
	..........
	..........
	..........
	ssbn.h ssbn.cpp
	)

* The neuron model uses burst length ('spb') as a variable. Add this variable to nest namespace.

	Add to	nestkernel/nest_names.h

	extern const Name spb;

	Add to nestkernel/nest_names.cpp

	const Name spb( "spb" );	
	


* Recompile and install nest
	make 
	make install

* Run the test script test_ssbn.py to check if module is installed. This also reproduces the Supp. Fig S1 (with FI curves)

* New module can also be installed using nestML, but this has not been tried yet.
More information about installing custom models in NEST:

https://nest.github.io/nest-simulator/extension_modules
https://github.com/nest/nestml




	
	

