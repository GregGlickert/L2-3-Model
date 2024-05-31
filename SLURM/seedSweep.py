import json


# something like this should work but I will let you do the rest
class EditJson:
    def __init__(self,json_file_path, key, new_parameter):
        self.json_file_path = json_file_path
        self.key = key
        self.new_parameter = new_parameter

    def EditJson(self):
    # Open the JSON file in read mode
        with open(self.json_file_path, 'r') as json_file:
            file_contents = json.load(json_file)
            # Update the value of the specified key
            if self.key in file_contents:
                file_contents[self.key] = self.new_parameter
            else:
                raise Exception(f"Key {self.key} not found in {self.json_file_path}")
                return
        # Write the updated JSON back to the file
        with open(self.json_file_path, 'w') as json_file:
            json.dump(file_contents, json_file, indent=4)
        
        print(f"JSON file '{self.json_file_path}' modified successfully.")

    # could add the save synaptic properties functions in this also those might be helpful here...
            
            