#
# Makefile
# dlilien, 2018-12-06 13:15
#
#

PYQT=pyuic5
UI_SOURCES=RawPickGUI.ui
UI_PY=$(UI_SOURCES:.ui=.py)

all: ui

ui: $(UI_PY)

$(UI_PY): %.py: %.ui
	$(PYQT)	-x $< -o $@
	sed 's/mplfigcanvaswidget/.mplfigcanvaswidget/' $@ | grep -v setShortcutVisibleInContextMenu > dum
	mv dum $@
# vim:ft=make
#
