from crystal_refinement.SHELXFile import SHELXFile

def main():
    test_file = "/Users/julialing/Documents/DataScience/crystal_refinement/4-2-1-4_single_crystal/Example_from_slides/7.ins"

    with open(test_file) as f:
        text = f.read()
        print text
        print "\n\n" + "~"*50 + "\n\n"
        file_obj = SHELXFile(text)

        print file_obj.write_ins()


if __name__ == "__main__":
    main()