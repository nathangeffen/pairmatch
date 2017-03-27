#/bin/sh
# Execute experiments for paper 2

echo Executing tests 1 to 5
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_01.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_02.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_03.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_04.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_05.csv &
wait

echo Executing tests 6 to 10
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_06.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_07.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_08.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_09.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_10.csv &
wait


echo Executing tests 11 to 15
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_11.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_12.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_13.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_14.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_15.csv &
wait
