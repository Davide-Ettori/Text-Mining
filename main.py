def read_data(data_file_path, parameters_file_path):
    T = list()
    MIS = dict()
    M = set()

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
                item = int(l[start:end])
                MIS[item] = float(l.split(' ')[-1])
    
    with open(data_file_path, 'r') as data_file:
        data_lines = data_file.readlines()
        for line in data_lines:
            temp = [int(x) for x in line.strip().split(',')]
            T.append(set(temp))
            for item in temp:
                M.add((item, get_MIS(item, MIS)))

    M = list(M)
    M.sort(key=lambda x: x[1])
    return T, MIS, SDC, [x[0] for x in M]

def write_output(output_file_path, F_k):
    with open(output_file_path, 'w') as output_file:
        for idx, list_itemset in enumerate(F_k, start=1):
            output_file.write(f"(Length-{idx} {len(list_itemset)}\n")

            for itemset in list_itemset:
                itemset_str = " ".join(str(item) for item in itemset.items)
                output_file.write(f"    ({itemset_str}) : {itemset.freq_count} : {itemset.tail_count}\n")
                #output_file.write(f"    ({itemset_str}) : {itemset.count}\n")
            output_file.write(")\n")

    return

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
        if supports[item] >= get_MIS(item, MIS):
            first_index = i
            L.append(item)
            base_support = supports[item] # the support that will be used for inserting the following items in L
            break
    
    if first_index == -1: # if no item satisfies the MIS, terminate
        print("No Item satisfies the MIS")
        return None, None
    
    for i in range(first_index + 1, len(M)): # insert the remaining items in L
        if supports[M[i]] >= base_support:
            L.append(M[i])
    
    return L, supports

def get_MIS(item, MIS):
    if item in MIS:
        return MIS[item]
    return MIS["rest"]

class Itemset():
    def __init__(self, items):
        self.items = items # must be sorted by MIS value
        self.items_set = set(items)
        self.freq_count = 0
        self.tail_count = 0
        
    def can_join(self, other_item, phi, supports): # return True iff the items can be joined
        return self.items[:-1] == other_item.items[:-1] and \
                    self.items[-1] < other_item.items[-1] and \
                        abs(supports[self.items[-1]] - supports[other_item.items[-1]]) <= phi
    
    def join(self, other_item): # return a new Itemset that is the join of the two  
        return Itemset(self.items + [other_item.items[-1]])
    
    def to_set(self):
        return self.items_set

def get_subsets(s):
    subsets = list()
    for i in range(len(s)):
        subsets.append(set(s[:i] + s[i + 1:]))
    return subsets

def lvl2_candidate_gen(L, SDC, MIS, supports):
    C2 = list()

    for i, l in enumerate(L):
        if supports[l] >= get_MIS(l, MIS):
            for h in L[i + 1:]:
                if supports[h] >= get_MIS(l, MIS) and abs(supports[h] - supports[l]) <= SDC:
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
                    if c.items[0] in s or get_MIS(c.items[0], MIS) == get_MIS(c.items[1], MIS):
                        if prune_candidates(s, F_k):
                            C_k.pop()   
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
    F_k.append([Itemset([item]) for item in L if supports[item] >= get_MIS(item, MIS)])

    # check freq vs tail count for F_1 itemsets
    for i in range(len(F_k[-1])): # adding the count attribute also for F_1 itemsets. From F_2 on is already done in the loop
        F_k[-1][i].freq_count = round(supports[F_k[-1][i].items[0]] * len(T))
        
    while (True):
        if k == 2:
            C_k = lvl2_candidate_gen(L, SDC, MIS, supports)
        else:
            C_k = MScandidate_gen(F_k[-1], SDC, MIS, supports)
        for t in T:
            for c in C_k:
                temp_set = c.to_set()
                if contains(temp_set, t):
                    c.freq_count += 1
                temp_set.remove(c.items[0])
                if contains(temp_set, t):
                    index = find_itemset(temp_set, F_k[-1])
                    F_k[-1][index].tail_count += 1
        F_k.append([c for c in C_k if c.freq_count/len(T) >= get_MIS(c.items[0], MIS)])
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
    output_file_path = 'out.txt'
    T, MIS, SDC, M = read_data(data_file_path, parameters_file_path)

    print(T)
    print(MIS)
    print(SDC)
    print(M)
    F_k = MSApriori(T, MIS, SDC, M)

    write_output(output_file_path, F_k)