from bioml import *
import os
import pandas as pd



def verify_bunch_models(folder_path):

    check_results = pd.DataFrame(columns=["Model Name", "Mass Action", "Reversible", "Plausible", "Error"])



    for file_name in os.listdir(folder_path):

        if not file_name.endswith('.xml'): continue    
        
        bioml = BioML()

        bioml.read_file(folder_path, file_name)

        mass_action: bool = None

        reversible: bool = None

        plausible: bool = None

        error: str = None

        try:


            if bioml.check_mass_action_kinetics(raise_error=True):

                mass_action = True

                if bioml.check_model_reversibility(raise_error=True):

                    reversible = True

                    if bioml.check_kinetic_constants_thermo_compatibility(raise_error=True):

                        plausible = True

                    else:

                        plausible = False

                else:

                    reversible = False

            else:

                mass_action = False

        except Exception as e:

            error = str(e)

        check_results.loc[len(check_results)] = [file_name, mass_action, reversible, plausible, error]

    excel_full_path = os.path.join(folder_path, "check_results.xlsx")
        
    check_results.to_excel(excel_full_path, index=False)

    save_message = f"\n{'*' * 30} The Excel file \"check_results.xlsx\" has been successfully saved to \"{folder_path}\" {'*' * 30}"

    utility.message_printer(save_message)