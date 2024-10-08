## ANSYS MODEL for building in Taufkirchen
This is the building model used for the MDSI project: [MSDI page](https://www.mdsi.tum.de/en/mdsi/forschung/foerderformate/anschubfinanzierung/seismicsafety/)

Author: Wei-Teng Kao (ge2gak@mytum.de)

### Run the model
1. Specify the path of folder and open the APDL 
2. Run the line: ```/input,'main','mac'``` or copy the script in main.mac directly then paste to the command window (optional).

### Files management
1. ANSYS_Building_model/PREP: Folder for preprossesor
    - [BuildPara_Var.mac](/PREP/BuildPara_Var.mac) : Material properties and parameter
    - [BuildGeo_TK.mac](/PREP/BuildGeo_TK.mac) : Geometry of the building
    - [SSI_LPM_para.mac](/PREP/SSI_LPM_para.mac) : Parameter of LPM for SSI 
    - [MAT27_elem.mac](/PREP/MAT27_elem.mac) : Assign the material 
    - properties to element for SSI 


2. ANSYS_Building_model/SOLU: Folder for solver
    - [Modal_Analysis.mac](/SOLU/Modal_Analysis.mac) : Modal analysis 


3. ANSYS_Building_model/POSTL Folder for postprossesor

### Parameter
- General:
  - ```bool_check_real_shape``` : Plot the building in mesh with defined thickness.
  
- Preprossesor:
  - ```bool_inner``` : Control whether build the inner wall
  - ```bool_stair``` : If both bool_inner = 1 and bool_stair = 1, then stair will be built
  - ```bool_SSI```   : Control soil-structure interection (MATRIX27)

- Solution: 
  - ```Solu_type```  : ANYTYPE (APDL) if 2 = Modal, 3 = Harmonic, 4 = Transient analysis




