def my_fun(param1):
    param1*=2
    return str(param1)

def main():
    write_out = open("out.txt", 'w')
    write_out.write("import sys\n")
    write_out.write("print >> sys.stderr, ('this is it')")
    write_out.close()

main()

