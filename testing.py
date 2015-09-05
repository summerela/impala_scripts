from xml.dom import minidom

xmldoc = minidom.parse('interpro.xml')
itemlist = xmldoc.getElementsByTagName('dbinfo')
print(len(itemlist))
print(itemlist[0].attributes['dbname'].value)
# for s in itemlist:
#     print(s.attributes['dbname'].value), (s.childNodes[0].nodeValue)


with open('test_out.csv', 'w') as fout:
      fieldnames = ['db', 'info']
      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
      writer.writeheader()
      for e in doc.xpath('dbinfo'):
         title, authors = e.xpath('author/text()'), e.xpath('title/text()')[0]
         writer.writerow({'title': titleValue, 'author': authors.join(';')})

#/opt/cloudera/parcels/CDH/lib/hadoop-mapreduce/hadoop-streaming.jar
