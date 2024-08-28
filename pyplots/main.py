import potential_energy_volume
import task2_energy 
import task3
import task4_sol
import task4_liq

def main():
    potential_energy_volume.plot_energy_pot_volume()
    task2_energy.plot_energy_kinetic()
    task3.plot_energy_task3()
    task3.plot_positions_task3()
    task3.plot_msd()
    task4_sol.plot_main_solid()
    task4_sol.plot_temps_pressures()
    task4_liq.plot_main_liq()
    task4_liq.plot_temps_pressures()

if __name__ == "__main__":
    main()
