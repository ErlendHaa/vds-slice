ACCOUNT_NAME=
SAS=

python defined.py defined.segy
SEGYImport --url "file://." --vdsfile defined_default.vds defined.segy
VDSCopy "defined_default.vds" "azureSAS://$ACCOUNT_NAME.blob.core.windows.net/testdata/defined/defined_default" --compression-method=None -d "Suffix=?$SAS"

# or import directly to the cloud (and see note about decompression)
