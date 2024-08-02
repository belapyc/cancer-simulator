import numpy as np
def align_to_diagnosis(df):
    # Group by patient and find the diagnosis point
    diagnosis_points = df[df['months_since_diagnosis'] == 0].groupby('sample_id')['time_entry_to_origin'].first()

    # Function to adjust time for each patient
    def adjust_time(group):
        patient_id = group['sample_id'].iloc[0]
        diagnosis_time = diagnosis_points.get(patient_id, np.inf)
        group['adjusted_time'] = group['time_entry_to_origin'] - diagnosis_time
        return group

    # Apply the adjustment to each patient's trajectory
    aligned_df = df.groupby('sample_id').apply(adjust_time).reset_index(drop=True)

    # Filter out pre-diagnosis data points and any rows with inf values
    aligned_df = aligned_df[(aligned_df['adjusted_time'] >= 0) & (aligned_df['adjusted_time'] != np.inf)]

    # Check that patient was diagnosed
    diagnosed_patients = aligned_df[aligned_df['months_since_diagnosis'] == 0]['sample_id'].unique()
    aligned_df = aligned_df[aligned_df['sample_id'].isin(diagnosed_patients)]

    return aligned_df