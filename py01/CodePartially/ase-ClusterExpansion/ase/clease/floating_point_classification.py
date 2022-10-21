
class FloatingPointClassifier(object):
    def __init__(self, num_decimals):
        self.num_dec = num_decimals
        self.lut = []

    def toJSON(self):
        json_fields = {
            "num_dec": self.num_dec,
            "lut": self.lut
        }
        return json_fields

    @staticmethod
    def fromJSON(json_fields):
        obj = FloatingPointClassifier(json_fields["num_dec"])
        obj.lut = json_fields["lut"]
        return obj

    def get(self, float_num):
        """Return the classification string of float_num."""
        tol = 2 * 10**(-self.num_dec)

        for entry in self.lut:
            if abs(entry["float"] - float_num) < tol:
                # Entry already exists in LUT
                return str(entry["classification"])

        # Create classification
        int_form = int(float_num * 10**self.num_dec)

        if self.num_dec >= 1:
            classification = str(int_form)[0] + 'p' + str(int_form)[1:]
        else:
            classification = str(int_form)
        new_entry = {"float": float_num,
                     "classification": classification}
        self.lut.append(new_entry)
        return new_entry["classification"]
