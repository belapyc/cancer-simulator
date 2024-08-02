import numpy as np
import pandas as pd


class LungCancerProgressionGenerator:
    def __init__(self, num_patients, years=15, doubling_time=134):
        self.num_patients = num_patients
        self.months = years * 12
        self.states = ['Healthy', 'Stage I', 'Stage II', 'Stage III', 'Stage IV', 'Death']
        self.doubling_time = doubling_time / 30  # Convert 134 days to months

        # Age-specific cancer occurrence rates (annual per 100,000) for ages 50-70
        self.cancer_rates = {
            (50, 54): 33.2, (55, 59): 79.9, (60, 64): 140.5,
            (65, 69): 198.4, (70, 74): 262.8
        }

        # Survival rates by stage
        self.survival_rates = {
            1: {'1year': 0.7112, '5year': 0.3533},
            2: {'1year': 0.4815, '5year': 0.2089},
            3: {'1year': 0.3495, '5year': 0.0632},
            4: {'1year': 0.1436, '5year': 0.0}
        }

        self.five_year_survival = {
            1: 0.65,
            2: 0.40,
            3: 0.15,
            4: 0.5
        }

    def get_cancer_rate(self, age):
        for (lower, upper), rate in self.cancer_rates.items():
            if lower <= age <= upper:
                return rate / 100000 / 12  # Convert to monthly rate
        return self.cancer_rates[(65, 69)] / 100000 / 12  # Default to 85+ rate

    def get_monthly_death_probability(self, stage, months_since_diagnosis):
        if stage == 0:  # Healthy
            return 0.0

        # if months_since_diagnosis <= 12:
        # yearly_survival = self.survival_rates[stage]['1year']
        # else:
        #     yearly_survival = self.survival_rates[stage]['5year']
        five_year_survival = self.five_year_survival[stage]

        monthly_survival = five_year_survival ** (1 / 60)
        return 1 - monthly_survival

    def diagnose_cancer(self, current_state, use_blood_test=False):
        if use_blood_test:
            return current_state  # Perfect detection with blood test
        else:
            # 70% chance to diagnose late stages if no blood test and 30% chance for early stages
            if current_state == 1 or current_state == 2:
                rand = np.random.random()
                if rand < 0.3:
                    return current_state
            elif current_state == 3 or current_state == 4:
                rand = np.random.random()
                if rand < 0.7:
                    return current_state
            return 0  # Cancer not detected

    def generate_patient_trajectory(self, patient_id, use_blood_test=False):
        trajectory = []
        current_state = 0  # Start in Healthy state
        age = 50  # Random starting age between 50 and 70
        cancer_time = 0  # Time since cancer onset
        months_since_diagnosis = -1  # -1 indicates cancer has not been diagnosed
        months_since_last_test = 0
        diagnosis_month = None
        doublings = 0  # Number of doublings needed for next stage

        for month in range(self.months):
            origin_state = current_state
            target_state = current_state  # Initialize target_state to current_state

            if current_state == 0:  # Healthy state

                cancer_prob = self.get_cancer_rate(age)

                if np.random.random() < cancer_prob:
                    target_state = 1  # Cancer starts at Stage I
                    cancer_time = 0
                    doublings = 0

                    diagnosed_state = self.diagnose_cancer(target_state)
                    if diagnosed_state > 0:
                        target_state = diagnosed_state
                        months_since_diagnosis = 0
                        diagnosis_month = month

            elif 1 <= current_state <= 4:  # Cancer stages I-IV
                if months_since_diagnosis >= 0:  # Only increment if cancer has been diagnosed
                    months_since_diagnosis += 1
                cancer_time += 1

                # Check for death
                death_prob = self.get_monthly_death_probability(current_state, months_since_diagnosis)
                if np.random.random() < death_prob:
                    target_state = 5  # Death state
                elif cancer_time >= self.doubling_time * (doublings + 1) and current_state < 4:
                    target_state = current_state + 1
                    cancer_time = 0  # Reset cancer time for the new stage
                    doublings = 0  # Reset doublings for the new stage
                    if months_since_diagnosis == -1:  # Only check for diagnosis if not yet diagnosed
                        diagnosed_state = self.diagnose_cancer(target_state, use_blood_test)
                        if diagnosed_state > 0:
                            # target_state = diagnosed_state
                            months_since_diagnosis = 0
                            diagnosis_month = month
                else:
                    doublings = int(cancer_time / self.doubling_time)  # Update doublings count

            # Blood test every 2 years (24 months)
            if use_blood_test and months_since_last_test >= 24 and months_since_diagnosis == -1:
                months_since_last_test = 0
                if target_state > 0:  # If cancer is present
                    months_since_diagnosis = 0
                    diagnosis_month = month

            if target_state != 0 or current_state != 0:  # Only record non-Healthy states
                trajectory.append({
                    'sample_id': patient_id,
                    'origin_state': origin_state,
                    'target_state': target_state,
                    'months_since_diagnosis': months_since_diagnosis,
                    'age_at_diagnosis': age - (month - diagnosis_month) / 12 if months_since_diagnosis >= 0 else None,
                    'time_entry_to_origin': month,
                    'time_transition_to_target': month + 1
                })

            if target_state == 5:  # If patient has died, stop tracking
                break

            current_state = target_state  # Update current_state for the next iteration
            age += 1 / 12  # Increase age by 1 month
            months_since_last_test += 1

        return trajectory

    def generate_dataset(self, use_blood_test=False):
        all_trajectories = []
        for i in range(self.num_patients):
            all_trajectories.extend(self.generate_patient_trajectory(i, use_blood_test))

        return pd.DataFrame(all_trajectories)