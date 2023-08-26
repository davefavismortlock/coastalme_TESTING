#!/bin/sh

cp in/test_suite/minimal_wave_angle_230/cme.ini .
./cme

cp in/test_suite/minimal_wave_angle_270/cme.ini .
./cme

cp in/test_suite/minimal_wave_angle_310/cme.ini .
./cme


cp in/test_suite/minimal_with_intervention_wave_angle_235/cme.ini .
./cme

cp in/test_suite/minimal_with_intervention_wave_angle_270/cme.ini .
./cme

cp in/test_suite/minimal_with_intervention_wave_angle_305/cme.ini .
./cme


cp in/test_suite/minimal_with_two_wavestations/cme.ini .
./cme

cp in/test_suite/minimal_with_sediment_input/cme.ini .
./cme


cp in/test_suite/Happisburgh/cme.ini .
./cme

#cp in/test_suite/Manuel_C003_0001/cme.ini .
#./cme

