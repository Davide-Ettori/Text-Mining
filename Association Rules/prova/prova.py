def read_data(data_file_path, parameters_file_path):
    T = list()
    MIS = dict()
    M = list()
    with open(data_file_path, 'r') as data_file:
        data_lines = data_file.readlines()
        for line in data_lines:
            temp = line.strip().split(',')
            T.append(set([int(x) for x in temp]))

    with open(parameters_file_path, 'r') as parameters_file:
        parameter_lines = parameters_file.readlines()
        for line in parameter_lines:
            l = line.strip()
            if "SDC" in l:
                SDC = float(l.split(' ')[-1])
            elif 'rest' in l:
                MIS["rest"] = float(l.split(' ')[-1])
            else:
                start = l.find('(') + 1
                end = l.find(')')
                item = l[start:end]
                MIS[int(item)] = float(l.split(' ')[-1])
                M.append((int(item), MIS[item]))
    
    M.sort(key=lambda x: x[1])
    return T, MIS, SDC, [x[0] for x in M]

def init_pass(M, T):
    supports = dict()
    
    L = list()
    for t in T: # count the number of occurences of each item in the transactions
        for item in t:
            if item in supports:
                supports[item] += 1
            else:
                supports[item] = 1

    for item in supports:
        supports[item] /= len(T) # convert the count to support (%)

    first_index = -1
    for i, item in enumerate(M): # find the first item that satisfies the MIS
        if supports[item] >= MIS[item]:
            first_index = i
            L.append(item)
            base_support = supports[item] # the support that will be used for inserting the following items in L
            break
    
    if first_index == -1: # if no item satisfies the MIS, terminate
        print("No Item satisfies the MIS")
        return None
    
    for i in range(first_index + 1, len(M)): # insert the remaining items in L
        if supports[M[i]] >= base_support:
            L.append(M[i])
    
    return L, supports

class Itemset():
    def __init__(self, items):
        self.items = items # must be sorted by MIS value
        self.count = 0
        
    def can_join(self, other_item, phi, supports): # return True iff the items can be joined
        return self.items[:-1] == other_item.items[:-1] and \
                    self.items[-1] < other_item.items[-1] and \
                        abs(supports[self.items[-1]] - supports[other_item.items[-1]]) <= phi
    
    def join(self, other_item): # return a new Itemset that is the join of the two  
        return Itemset(self.items + [other_item.items[-1]])
    
    def to_set(self):
        return set(self.items)

def get_subsets(s):
    subsets = list()
    for i in range(len(s)):
        subsets.append(set(s[:i] + s[i + 1:]))
    return subsets

def lvl2_candidate_gen(L, SDC, MIS, supports):
    C2 = list()

    for i, l in enumerate(L):
        if supports[l] >= MIS[l]:
            for h in L[i + 1:]:
                if supports[h] >= MIS[l] and abs(supports[h] - supports[l]) <= SDC:
                    C2.append(Itemset([l, h]))
    return C2

def prune_candidates(s, F_k):
    for f in F_k:
        if f.to_set() == s:
            return False
    return True

def MScandidate_gen(F_k, SDC, MIS, supports):
    C_k = []
    for f1 in F_k:
        for f2 in F_k:
            if f1.can_join(f2, SDC, supports):
                c = f1.join(f2)
                C_k.append(c)
                subsets = get_subsets(c.items)
                for s in subsets:
                    if s.contains(c.items[0]) or MIS[c.items[0]] == MIS[c.items[1]]:
                        if prune_candidates(s, F_k):
                            C_k.pop(c)   
                            break  
    return C_k    

def find_itemset(itemset, F_k):
    for i, f in enumerate(F_k):
        if f.to_set() == itemset:
            return i
    return -1

def MSApriori(T, MIS, SDC, M):
    L, supports = init_pass(M, T) # returns both the L list and the supports for each singular item
    k = 1
    F_k = list() 
    F_k.append([Itemset([item]) for item in L if supports[item] >= MIS[item]])
    while (True):
        if k == 2:
            C_k = lvl2_candidate_gen(L, SDC, MIS, supports)
        else:
            C_k = MScandidate_gen(F_k[-1], SDC, MIS, supports)
        for t in T:
            for c in C_k:
                if contains(c.to_set(), t):
                    c.count += 1
                if contains(c.to_set().remove(c.items[0]), t):
                    index = find_itemset(c.to_set().remove(c.items[0]), F_k[-1])
                    F_k[-1][index].count += 1
        F_k.append([c for c in C_k if c.count/len(T) >= MIS[c.items[0]]])
        if not F_k[-1]:
            break
        k += 1
    F_k.pop()
    return F_k


def contains(a, b):
    return a.issubset(b)

if __name__ == "__main__":
    data_file_path = 'data.txt'
    parameters_file_path = 'parameters.txt'
    T, MIS, SDC, M = read_data(data_file_path, parameters_file_path)

    print(T)
    print(MIS)
    print(SDC)
    print(M)
    MSApriori(T, MIS, SDC, M)