%% State Table
State_vector = {'x'; 'y'; 'alpha'; 'varphi_r'; 'varphi_l'; 'varphi_p';'x_dot';'y_dot';'alpha_dot';'varphi_dot_r';'varphi_dot_l';'varphi_dot_p'};
State_experiment4_t0 = [states(1,1), states(1,2), states(1,3), states(1,4), states(1,5), states(1,6), states(1,7), states(1,8), states(1,9), states(1,10), states(1,11), states(1,12)]';
State_experiment4_t40 = [states(40001,1), states(40001,2), states(40001,3), states(40001,4), states(40001,5), states(40001,6), states(40001,7), states(40001,8), states(40001,9), states(40001,10), states(40001,11), states(40001,12)]';

Ttable_states = table(State_vector, State_experiment4_t0, State_experiment4_t40)
%% Energy Table
Energy = {'Total Kinetic Energy'; 'Rotational Kinetic Energy of Otbot'; 'Translational Kinetic Energy of Otbot'; 'Rotational Kinetic Energy of the Chassis'; 'Translational Kinetic Energy of the Chassis'; 'Chassis Body (No wheels) Rotational Energy'; 'Right Wheel Rotational Energy'; 'Left wheel Rotational Energy'; 'Platform Rotational Energy'; 'Chassis Body (No Wheels) Translational Energy'; 'Right Wheel Translational Energy'; 'Left Wheel Translational Energy'; 'Platform Translational Energy'};
Energies_experiment4_t0 = [kinenvalues(1),kinenrot_sys(1),kinentra_sys(1),kinenrot_chas(1),kinentra_chas(1),kinenrot_b(1),kinenrot_r(1),kinenrot_l(1),kinenrot_p(1),kinentra_b(1),kinentra_r(1),kinentra_l(1),kinentra_p(1)]';
Energies_experiment4_t40 = [kinenvalues(40001),kinenrot_sys(40001),kinentra_sys(40001),kinenrot_chas(40001),kinentra_chas(40001),kinenrot_b(40001),kinenrot_r(40001),kinenrot_l(40001),kinenrot_p(40001),kinentra_b(40001),kinentra_r(40001),kinentra_l(40001),kinentra_p(40001)]';

Ttable_Energies = table(Energy, Energies_experiment4_t0, Energies_experiment4_t40)

%% Save Tables
writetable(Ttable_states,'Table_of_States.xlsx')
writetable(Ttable_Energies,'Table_of_Energies.xlsx')

