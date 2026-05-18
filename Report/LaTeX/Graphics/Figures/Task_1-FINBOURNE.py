#%%

# %%

from collections import Counter
import numpy as np

vacations = [1,2,1,3,4,1,2,5,1,2,3,4,1,2]

n_locations = list(Counter(vacations).keys()) # equals to list(set(words))
values = list(Counter(vacations).values()) # counts the elements' frequency

n_check_array = len(n_locations)

temp_check_vacations = vacations[0:n_check_array]

while all(elem in temp_check_vacations for elem in n_locations) == False:
    print(n_check_array)
    for i in range(len(vacations)-n_check_array+1):
        temp_check_vacations = vacations[i:i+n_check_array]
        print(temp_check_vacations)
        # if all(elem in temp_check_vacations for elem in n_locations):
        # print(len(temp_check_vacations))
    n_check_array += 1



#%%----------------------------------------------------------------

from collections import Counter
import numpy as np

vacations = [2,1,1,3,2,1,1,3,4,5,3,2,1,6,6,6,6,9]

"""
def Solution(A):

    n_locations = list(Counter(vacations).keys()) # equals to list(set(words))
    values = list(Counter(vacations).values()) # counts the elements' frequency

    n_check_array = len(n_locations)

    for i in range(len(vacations)-n_check_array+1):
        temp_check_vacations = vacations[i:i+n_check_array]
        if all(elem in temp_check_vacations for elem in n_locations):
            return len(temp_check_vacations)
        if i == max(range(len(vacations)-n_check_array)):
            n_check_array += 1

Solution(vacations)
"""

"""
n_locations = list(Counter(vacations).keys()) # equals to list(set(words))
values = list(Counter(vacations).values()) # counts the elements' frequency

n_check_array = len(n_locations)

temp_check_vacations = vacations[i:i+n_check_array]

while all(elem in temp_check_vacations for elem in n_locations) == False:
    for i in range(len(vacations)-n_check_array+1):
        temp_check_vacations = vacations[i:i+n_check_array]
    if all(elem in temp_check_vacations for elem in n_locations):
        print(len(temp_check_vacations))
    if i == max(range(len(vacations)-n_check_array)):
        n_check_array += 1

"""
#%%

from collections import Counter
import numpy as np

vacations = [1,2,1,3,4,1,2,5,1,2,3,4,1,2,7]

n_locations = list(Counter(vacations).keys()) # equals to list(set(words))
values = list(Counter(vacations).values()) # counts the elements' frequency

n_check_array = len(n_locations)

temp_check_vacations = vacations[0:n_check_array]

while all(elem in temp_check_vacations for elem in n_locations) == False:
    print(n_check_array)
    for j in range(len(vacations)-n_check_array+1):
        temp_check_vacations = vacations[j:j+n_check_array]
        print(temp_check_vacations)
        # if all(elem in temp_check_vacations for elem in n_locations):
        # print(len(temp_check_vacations))
    n_check_array += 1
pass


# %%
