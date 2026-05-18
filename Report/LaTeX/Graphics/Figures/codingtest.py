#%%

# you can write to stdout for debugging purposes, e.g.
# print("this is a debug message")

N = 480185961

print(solution(N))

def solution(N):
    # write your code in Python 3.6

    a = str(bin(N))
    c = '1'
    gaps = []
    num_space = []

    pos_ones = [pos for pos, char in enumerate(a) if char == c]

    for i in range(len(pos_ones)-1):
        if pos_ones[i+1] == pos_ones[i] + 1:
            pass
        elif pos_ones[i+1] != pos_ones[i] + 1:
            gaps.append([pos_ones[i], pos_ones[i+1]])

    if len(gaps) == 0:
        print('There are no binary gaps')
        return 0

    for i in range(len(gaps)):
        space = gaps[i]
        num_space.append(space[1] - space [0] - 1)

    longest_gap = max(num_space)
    
    print('Number of binary gaps are', len(gaps))
    print('The respective gap sizes are', num_space)


    return longest_gap





                

# %%

# %%
