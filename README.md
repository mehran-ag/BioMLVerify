# BioMLVerify
This code is a tool to verify biological models written in XML format.

The code investigates mass conservation, charge conservation and thermodynamic compatibility of reaction networks in SBML and CellML models.

Currently, the code can check reactions governed by Mass Action law.

If the species in the reaction are annotated by ChEBI codes, mass and charge conservations can be checked automatically. Otherwise, there is an alternative method os manually assigning chemical elements or chemical moeities for mass conservation check. Charge can also be added to charged species if they have a charge, otherwise it will be considered zero.

To use the tool, an instance of BioML object must be created.

Using BioML instance, the model can be read and verified easily by only using two functions: read_file and verify_model.
