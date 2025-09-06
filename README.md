# BioMLVerify

**BioML: A Tool for Verifying Biological Models**

This tool verifies biological models written in XML format, focusing on:

- **Mass conservation**  
- **Charge conservation**  
- **Thermodynamic compatibility** of reaction networks  

It supports models in both **SBML** and **CellML** formats.

---

## Supported Reaction Types
- Currently, the tool can check reactions governed by the **Law of Mass Action**.

---

## Conservation Checks
- **Automatic check**:  
  If reaction species are annotated with **ChEBI codes**, mass and charge conservation are checked automatically.  

- **Manual check**:

  If ChEBI codes are not available, chemical elements or moieties can be assigned manually for mass conservation.  
  Charges can also be assigned manually; unassigned species are assumed to have zero charge.  

---

## Usage

To use the tool:

1. Create an instance of the `BioML` object:

   ```python
   bioml = BioML()

2. Use just two functions to read and verify a model:

```python
bioml.read_file(folder_path="C:/my folder", file_name="model_file.xml")
# or
bioml.read_file(folder_path="C:/my folder", file_name="model_file.cellml")

bioml.verify_model(printing=True)