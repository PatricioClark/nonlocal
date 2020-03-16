LIB_DIR = lib

default: spherical

spherical: setup.py spherical.pyx #$(LIB_DIR)/libspherical.a
	python3 setup.py build_ext --inplace && rm -f spherical.c && rm -Rf build

# $(LIB_DIR)/libspherical.a:
# 	make -C $(LIB_DIR) libspherical.a

clean:
	rm *.so
