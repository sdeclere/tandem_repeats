from PIL import Image as Img
from math import ceil
import sys 
import os 

def to_pdf(aline):
    png_files = ['cache/%s.png' % x for x in aline[1].split(',') ]
    png_files = list(filter(os.path.isfile, png_files))
    print( aline )
    try:
        imgs = [Img.open(i) for i in png_files]
        (width, height) = imgs[0].size
        #print((width, height))
        #(width, height) = 595 , 842
        nb_lines = ceil(len(imgs)/3)
        print (len(imgs))
        print (nb_lines)
        new_img = Img.new('RGB', (3 * width, nb_lines * height), color="white")
        i=0
        for j in range(nb_lines):
            for k in range(3):
                if i >= len(imgs):
                    pass 
                else: 
                    new_img.paste(im=imgs[i], box=(k * width, j * height))
                    i+=1

        #new_img.save(os.path.join('book', 'toto' + '.pdf'), format='pdf')
    except Exception as e:
        print(e)
    return new_img

with open (sys.argv[1]) as fd:
    for line in fd:
        sl = line.strip().split('|')

        if sl[0] == sys.argv[2]:
        #if sl[5] != 'na':
            output_filename = 'book/%s.pdf' % (sl[0])
            print ('generating %s ...' % output_filename)
            fam_members=sl[1].split(',')
            if len(fam_members) > 100:
                pages =[]
                for i in range(0, len(fam_members), 100):
                    page = sl[:]
                    page[1]= ','.join(fam_members[i:i+100])
                    pages.append( to_pdf(page) )
                first = pages.pop(0)
                first.save(output_filename, save_all=True, append_images=pages)
            else:
                img = to_pdf(sl)
                img.save(output_filename, format='pdf')

