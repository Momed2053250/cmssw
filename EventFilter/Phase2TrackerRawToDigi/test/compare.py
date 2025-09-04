import re

def parse_file(filename):
    data = {
        'header': [],
        'offset': [], 
        'lines': [],
        'strip': [],
        'pixel': [],
        'is2SModule': [],
        'channelOffset16': [],
        'idx': [],
        'numStripClusters': [],
        'numPixelClusters': []
    }
    
    current_section = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            if 'headerWords' in line and '=' in line:
                try:
                    value = int(line.split('=')[1].strip())
                    data['header'].append(value)
                except:
                    pass
                    
            elif 'offsetWords' in line and '=' in line:
                try:
                    value = int(line.split('=')[1].strip())
                    data['offset'].append(value)
                except:
                    pass
                    
            elif 'Lines[' in line and '=' in line:
                try:
                    value = int(line.split('=')[1].strip())
                    data['lines'].append(value)
                except:
                    pass
                    
            elif 'Strip Cluster words:' in line:
                current_section = 'strip'
            elif 'Pixel Cluster words:' in line:
                current_section = 'pixel'
                
            elif current_section and line.replace(' ', '').isdigit():
                numbers = [int(x) for x in line.split()]
                data[current_section].extend(numbers)
            
            # Parse the additional printf statements
            elif 'is2SModule is:' in line:
                try:
                    value = int(line.split(':')[1].strip())
                    data['is2SModule'].append(value)
                except:
                    pass
                    
            elif 'ChannelOffset16 is:' in line:
                try:
                    value = int(line.split(':')[1].strip())
                    data['channelOffset16'].append(value)
                except:
                    pass
                    
            elif 'idx is:' in line:
                try:
                    value = int(line.split(':')[1].strip())
                    data['idx'].append(value)
                except:
                    pass
                    
            elif 'n strip clusters are:' in line:
                try:
                    value = int(line.split(':')[1].strip())
                    data['numStripClusters'].append(value)
                except:
                    pass
                    
            elif 'n pixel clusters are:' in line:
                try:
                    value = int(line.split(':')[1].strip())
                    data['numPixelClusters'].append(value)
                except:
                    pass
                
    return data

# Parse both files
gpu = parse_file('GPUcheck.log')
serial = parse_file('SerialOnlycheck.log')

print("="*60)
print("CLUSTER WORD COUNTS")
print("="*60)
print(f"HeaderWords count - GPU: {len(gpu['header'])}, Serial: {len(serial['header'])}")
print(f"OffsetWords count - GPU: {len(gpu['offset'])}, Serial: {len(serial['offset'])}")
print(f"Lines count - GPU: {len(gpu['lines'])}, Serial: {len(serial['lines'])}")
print(f"Strip clusters - GPU: {len(gpu['strip'])}, Serial: {len(serial['strip'])}")
print(f"Pixel clusters - GPU: {len(gpu['pixel'])}, Serial: {len(serial['pixel'])}")

print("\n" + "="*60)
print("ADDITIONAL DEBUG COUNTS")
print("="*60)
print(f"is2SModule count - GPU: {len(gpu['is2SModule'])}, Serial: {len(serial['is2SModule'])}")
print(f"channelOffset16 count - GPU: {len(gpu['channelOffset16'])}, Serial: {len(serial['channelOffset16'])}")
print(f"idx count - GPU: {len(gpu['idx'])}, Serial: {len(serial['idx'])}")
print(f"numStripClusters count - GPU: {len(gpu['numStripClusters'])}, Serial: {len(serial['numStripClusters'])}")
print(f"numPixelClusters count - GPU: {len(gpu['numPixelClusters'])}, Serial: {len(serial['numPixelClusters'])}")
