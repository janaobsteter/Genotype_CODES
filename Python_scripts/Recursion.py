# -*- coding: utf-8 -*-
#globina(oseba) - kako globoka je njgeova rodbina - če nima otrk = 1
#razmišljaj: 1 + globina mojih otrok, poglej le zase, če imaš otroke
#zadnji ma globino 1

def globina(oseba):
    potomci = otroci(oseba)
    if not potomci:
        glob = 1
    else:
        glob =+ 1
        for otrok in potomci:
            globina
    return glob 
        

# -*- coding: utf-8 -*-
rodbina = {
    "Adam": ["Matjaz", "Cilka", "Daniel"],
    "Aleksander": [],
    "Alenka": [],
    "Barbara": [],
    "Cilka": [],
    "Daniel": ["Elizabeta", "Hans"],
    "Erik": [],
    "Elizabeta": ["Ludvik", "Jurij", "Barbara"],
    "Franc": [],
    "Herman": ["Margareta"],
    "Hans": ["Herman", "Erik"],
    "Jozef": ["Alenka", "Aleksander", "Petra"],
    "Jurij": ["Franc", "Jozef"],
    "Ludvik": [],
    "Margareta": [],
    "Matjaz": ["Viljem"],
    "Petra": [],
    "Tadeja": [],
    "Viljem": ["Tadeja"],
}

def otroci(oseba):
    return rodbina[oseba]


import unittest
class Test(unittest.TestCase):
    def test_globina(self):
        self.assertEqual(globina("Alenka"), 1)
        self.assertEqual(globina("Cilka"), 1)
        self.assertEqual(globina("Viljem"), 2)
        self.assertEqual(globina("Jurij"), 3)
        self.assertEqual(globina("Hans"), 3)
        self.assertEqual(globina("Elizabeta"), 4)
        self.assertEqual(globina("Adam"), 6)

    def test_globina_imena(self):
        self.assertEqual(globina_imena("Alenka", "Alenka"), 0)
        self.assertIsNone(globina_imena("Matjaz", "Erik"))
        self.assertEqual(globina_imena("Matjaz", "Tadeja"), 2)
        self.assertEqual(globina_imena("Matjaz", "Viljem"), 1)
        self.assertEqual(globina_imena("Matjaz", "Matjaz"), 0)
        self.assertEqual(globina_imena("Matjaz", "Viljem"), 1)
        self.assertEqual(globina_imena("Daniel", "Ludvik"), 2)
        self.assertEqual(globina_imena("Daniel", "Alenka"), 4)
        self.assertEqual(globina_imena("Hans", "Erik"), 1)

    def test_globje_od_n(self):
        self.assertEqual(globje_od_n("Jozef", 0), {"Jozef", "Alenka", "Petra", "Aleksander"})
        self.assertEqual(globje_od_n("Jozef", 1), {"Alenka", "Petra", "Aleksander"})
        self.assertEqual(globje_od_n("Jozef", 2), set())

        self.assertEqual(globje_od_n("Adam", 0), set(rodbina))
        self.assertEqual(globje_od_n("Adam", 1), set(rodbina) - {"Adam"})
        self.assertEqual(globje_od_n("Adam", 2), set(rodbina) - {"Adam"} - set(rodbina["Adam"]))
        self.assertEqual(globje_od_n("Adam", 4), {"Franc", "Jozef", "Margareta", "Alenka", "Petra", "Aleksander"})


if __name__ == "__main__":
    unittest.main()
