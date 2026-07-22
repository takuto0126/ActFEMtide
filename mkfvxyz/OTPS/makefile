all:predict_tide extract_HC extract_local_model

predict_tide: predict_tide.f90 subs.f90 constit.h
	gfortran -o predict_tide -fconvert=swap -frecord-marker=4 predict_tide.f90 subs.f90
extract_HC:  extract_HC.f90 subs.f90
	gfortran -o extract_HC -fconvert=swap -frecord-marker=4 extract_HC.f90 subs.f90
extract_local_model: extract_local_model.f90 subs.f90
	gfortran -o extract_local_model -fconvert=swap -frecord-marker=4 extract_local_model.f90 subs.f90

clean:
	rm -f predict_tide
	rm -f extract_HC
	rm -f extract_local_model


