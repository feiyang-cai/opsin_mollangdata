[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

OPSIN - Open Parser for Systematic IUPAC Nomenclature
=====================================================

__License: [MIT License](https://opensource.org/licenses/MIT)__  

This repository is a customized OPSIN fork used for anonymous review. 

This repository is provided primarily as **code reference**. For the end-to-end MolLangData data-generation workflow and instructions on how this fork is used, please refer to the companion MolLangData link:

<https://anonymous.4open.science/r/MolLangData>

OPSIN is a Java library for IUPAC name-to-structure conversion supporting SMILES, CML, and InChI-related workflows.

Java 8 or higher is required.

### Simple Usage Examples

#### Convert a chemical name to SMILES

`java -jar opsin-cli-2.8.0-jar-with-dependencies.jar -osmi input.txt output.txt`

where `input.txt` contains one chemical name per line.

```java
NameToStructure nts = NameToStructure.getInstance();
String smiles = nts.parseToSmiles("acetamide");
```

### Availability

This project can be built from source and produces:

- `opsin-cli-<version>-jar-with-dependencies.jar`
- `opsin-core-<version>-jar-with-dependencies.jar`
- `opsin-inchi-<version>-jar-with-dependencies.jar`

### Building from Source

Run from the repository root:

```bash
mvn package
```

Artifacts are produced under:

- `opsin-cli/target`
- `opsin-core/target`
- `opsin-inchi/target`

### Customized Fork Notes

This fork includes custom functionality used by the companion dataset-generation workflow.

- This repository is intended mainly for code reference.
- For dataset-generation usage and workflow details, refer to the MolLangData link above.

### About OPSIN

OPSIN supports a broad range of systematic nomenclature, including:

- alkanes, alkenes, alkynes, and heteroatom chains,
- IUPAC 1993 recommended rings,
- Hantzsch-Widman systems,
- spiro and fused ring systems,
- bridge prefixes,
- stereochemical descriptors,
- functional replacement nomenclature,
- many functional classes and biomolecule-related naming patterns.

Some less common stereochemical terms and natural-product-specific operations remain unsupported.

### License

This project is licensed under the MIT License.
