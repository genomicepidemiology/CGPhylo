# Define the function to convert storage sizes to bytes
def convert_to_bytes(size_str):
    # Define unit multipliers
    unit_multipliers = {
        'K': 1024,
        'M': 1024 ** 2,
        'G': 1024 ** 3,
        'T': 1024 ** 4,
    }

    # Extract the number and unit from the string
    num, unit = size_str[:-1], size_str[-1]
    return float(num) * unit_multipliers.get(unit, 1)

# List of storage sizes
storage_sizes = [
    "272G", "32K", "32K", "72K", "32K", "9.8G", "32K", "46G", "55G", "619G", "32K",
    "176M", "32K", "514G", "32K", "32K", "32K", "23G", "32K", "32K", "755G", "32K",
    "32K", "32K", "32K", "32K", "32K", "32K", "32K", "32K", "2.9G", "32K", "278G",
    "975G", "222G", "644G", "32K", "32K", "2.5T", "466G", "3.7T", "32K", "346G", "70G",
    "33G", "32K", "32K", "32K", "32K", "32K", "32K", "32K", "91G", "234G", "841G", "32K",
    "32K", "1.1T", "32K", "5.9T", "32K", "32K", "11M", "64K", "32K", "32K", "32K", "32K",
    "13G", "25G", "1.4G", "132G", "64K", "2.7M", "69G", "32K", "32K", "32K", "32K", "455M",
    "136K", "163G", "289G", "32K", "32K", "2.8M", "32K", "5.3G", "32K", "32K", "7.1G", "936G",
    "5.5G", "32K", "32K", "806G", "32K", "32K", "118G", "32G", "32K", "32K", "59G", "32K",
    "346G", "32K", "11G", "32K", "64G", "3.3T", "32K", "421G", "32K", "32K", "375M", "32K",
    "32K", "32K", "100M", "32K", "104K", "81G", "2.5G", "104K", "32K", "1.8T", "852G", "32K",
    "944G", "32K", "32K", "32K", "2.4T", "242G", "256G", "2.3T", "32K", "32K", "32K", "1.3T",
    "32K", "512G", "32K", "32K", "4.6G", "32K", "3.0T", "75G", "95G", "165G", "32K", "624G",
    "535G", "40G", "930G", "32K", "280K", "181G", "176M", "32K", "512M", "46G", "146G", "147G",
    "424K", "520K", "32K", "32K", "489G", "32K", "198G", "222M", "32K", "32K", "32K", "26G",
    "31G", "32K", "4.2G", "32K", "88G", "221M", "32K", "32K", "32K", "32K", "104K", "32K",
    "32K", "32K", "11M", "642G", "106G", "32K", "78G", "32K", "32K", "4.6T", "2.3T", "32K",
    "32K"
]

# Correct the syntax error and complete the calculation
total_bytes = sum(convert_to_bytes(size) for size in storage_sizes)

# Convert the total bytes back to a more readable format
def bytes_to_readable(size_bytes):
    if size_bytes < 1024:
        return f"{size_bytes} bytes"
    elif size_bytes < 1024 ** 2:
        return f"{size_bytes / 1024:.2f} KB"
    elif size_bytes < 1024 ** 3:
        return f"{size_bytes / (1024 ** 2):.2f} MB"
    elif size_bytes < 1024 ** 4:
        return f"{size_bytes / (1024 ** 3):.2f} GB"
    else:
        return f"{size_bytes / (1024 ** 4):.2f} TB"

# Convert total bytes to a readable format
total_readable = bytes_to_readable(total_bytes)
print (total_readable)
