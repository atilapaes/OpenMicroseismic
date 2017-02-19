## Synopsis:

This project aims to organize ideas of the OpenMicroseismic project workflow.<br/>
Suggestions are welcome in its field. The developers are supposed to review them periodically and incorporate in the proper place.<br/>
Developers of subsequent packs must specify inputs and outputs for each pack and match the borders.<br/>

## Software structure and file organization:

By default, projects in Python use the convention “Snake Case”. In short, the variables and files are named in lowercase and words are separated by underscore.<br/>
In order to organize the project, the following structure is adopted: “project_pack_method_file_version”

For each pack, there will be 3 files:<br/>
Main: the main project. Composed the by the outer processing methods like general organization and multiprocessing tools.<br/>
Functions: All the functions used in the pack. Preferentially with comprehensive names, well documented and clear definition of Input and Output formats.<br/>
Parameters: Definition of all parameters used in the pack. Preferentially blocked in categories.<br/>

Example: <br/>
Project: OpenMicroseismic<br/>
Pack: Potential-Event Detection<br/>
Method: Energy-stack<br/>
The Potential-event detection pack (Version 2) is composed by the following files:<br/>
om_ped_es_main_v2.py<br/>
om_ped_es_functions_v2.py<br/>
om_ped_es_parameters_v2.py<br/>

Exceptions:<br/>
Functions with useful capabilities can be stored separately by the name:<br/>
om_general_function_name.py 

## Motivation

"Together we can achieve more"

## Contributors

Before you start uploading files to this folder, please add another folder with your name in and upload your files inside this folder with your name.

## License
