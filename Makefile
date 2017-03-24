pkg = $(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)
dir = '../$(pkg)'

document:
	R -e "devtools::document($(dir))"

check:
	R -e "devtools::check($(dir), document = FALSE)"

install: document
	R CMD INSTALL --no-multiarch --with-keep.source $(dir)

cleanrebuild: document
	R CMD INSTALL --preclean --no-multiarch --with-keep.source $(dir)
