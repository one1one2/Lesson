MAKEFILE = makefile

default: doc program 
	@cd doc && $(MAKE)
	@cd program && $(MAKE)

all: doc program 
	@cd doc && $(MAKE) all
	@cd program && $(MAKE) all
	@cd images && $(MAKE) all

clean:
	@cd doc && $(MAKE) clean
	@cd program && $(MAKE) clean
