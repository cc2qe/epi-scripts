#!/usr/bin/env python

import sys

def filter_impact(values):
    impact_order = ["synonymous_variant", "missense_variant", "stop_gained", "stop_lost", "start_lost"]
    
    # Filter out unwanted values
    filtered_values = [value for value in values if value in impact_order]
    
    # Sort the filtered values based on impact_order and get the highest impact
    sorted_values = sorted(filtered_values, key=lambda x: impact_order.index(x))
    highest_impact = sorted_values[0] if sorted_values else None
    
    return highest_impact

def main():
    # Print the header
    print("CHR_A\tBP_A\tSNP_A\tMAF_A\tCHR_B\tBP_B\tSNP_B\tMAF_B\tR\tD\tGENE_A\tGENE_B\tIMPACT_A\tIMPACT_B\tCADD_A\tCADD_B")

    for line in sys.stdin:
        # Split the tab-delimited line into columns
        columns = line.strip().split('\t')
        
        # Check if there are enough columns
        if len(columns) >= 14:
            # Extract values from columns 13 and 14
            impact_a_values = columns[12].split('&')
            impact_b_values = columns[13].split('&')
            
            # Filter and get the highest impact for each column
            filtered_impact_a = filter_impact(impact_a_values)
            filtered_impact_b = filter_impact(impact_b_values)
            
            # Check if both values are not None before printing the line
            if filtered_impact_a is not None and filtered_impact_b is not None:
                # Update columns with filtered values
                columns[12] = filtered_impact_a
                columns[13] = filtered_impact_b
                
                # Print the modified line
                print('\t'.join(columns))

if __name__ == "__main__":
    main()
