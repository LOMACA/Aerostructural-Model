
class AirfoilData:
    def __init__(self):
        pass

    def data_extraction(self, file_path, angle_start, angle_end):
        data = []
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines[12:]:  
                    if not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) == 7:
                        alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr = map(float, parts)
                        if angle_start <= alpha <= angle_end:
                            data.append({'alpha': alpha, 'CL': CL, 'CD': CD, 'CDp': CDp, 'CM': CM, 'Top_Xtr': Top_Xtr,
                                         'Bot_Xtr': Bot_Xtr})
            return data
        except Exception as e:
            print(f"Failed to read the file: {e}")
            return None
