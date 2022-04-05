#!/usr/bin/env python

import obja
import numpy as np
import sys
from os.path import isfile

class Decimater(obja.Model):
    def __init__(self):
        super().__init__()
        self.batch_vertices = []
        self.batch_faces = []
        self.K = {}

    def contract(self, output):
        
        operations = []
        n = 100 # nombre d'edge a collapse a chaque iteration
        # vertex = [coordonnees] (indice implicite par la liste)
        # face = [indice1, indice2, indice3] (indice implicite par la liste)

        self.batch_vertices = [(ind,vert) for ind,vert in enumerate(self.vertices)]
        self.batch_faces = [(ind,face) for ind,face in enumerate(self.faces)]
        self.K = self.compute_K()

        last_nb_vertices = 0

        while (len(self.batch_vertices) != last_nb_vertices):
            # Store indices to get the correspondance between batch ind and global ind
            vert_ind_2_batch_vert_ind = {i:batch_ind for batch_ind, (i,v) in enumerate(self.batch_vertices)}

            print(len(self.batch_vertices))

            # edges_with_error = [(edge, error), ...], edge = (v1,v2)
            edges_with_error = self.compute_error() # retourne les erreurs associées aux edges, et les edges
            potential_edge_collapse = dict(sorted(edges_with_error.items(), key=lambda item: item[1])).keys() # ordonner les edges par erreur
            edges_to_collapse = self.select_valid_edges(potential_edge_collapse,n) 

            removed_vert_indices = []

            # on calcule les operations necessaires pour obj : suppression de v1
            for (v1,v2) in edges_to_collapse:
                removed_faces_indices = []
                # operations.append(('##', "Suppression de " + str(v1) + " (" + str(self.vertices[v1]) + ") et " + str(v2), ""))

                for batch_face_ind, (face_index, face) in enumerate(self.batch_faces):
                    if v1 in [face.a,face.b,face.c]:
                        # si v1 et v2 dans la face -> supprimer la face
                        if v2 in [face.a,face.b,face.c]:
                            removed_faces_indices.append(batch_face_ind)
                            # operations.append(('#', str(v1) + " et " + str(v2) + " meme face", ""))
                            operations.append(('f', face_index, face))
                        # si uniquement v1 dans la face, remplacer v1 par v2
                        else:
                            # operations.append(('#', "Uniquement v1 (" + str(v1) + ")", ""))
                            operations.append(('ef', face_index, obja.Face.clone(face)))

                            if v1 == face.a:
                                face.a = v2
                            elif v1 == face.b:
                                face.b = v2
                            else:
                               face.c = v2

                # Delete the vertex
                operations.append(('v', v1, self.vertices[v1]))
                removed_vert_indices.append(vert_ind_2_batch_vert_ind[v1])

                # Remove the faces
                for face_ind in sorted(removed_faces_indices, reverse=True):
                    self.batch_faces.pop(face_ind)
            
            last_nb_vertices = len(self.batch_vertices)

            # Removes the vertices
            for vert_ind in sorted(removed_vert_indices, reverse=True):
                self.batch_vertices.pop(vert_ind)

        deleted_faces = set()

        # Decimate the remaining vertices
        for (vertex_index, vertex) in self.batch_vertices:
            for (face_index, face) in self.batch_faces:
                if face_index not in deleted_faces:
                    if vertex_index in [face.a,face.b,face.c]:
                        deleted_faces.add(face_index)
                        # operations.append(('#', "fin", ""))
                        operations.append(('f', face_index, obja.Face.clone(face)))
            operations.append(('v', vertex_index, vertex))


        # print("==== Ecriture dans le fichier ====")

        # for (ty, index, value) in operations:
        #     if ty[0] == '#':
        #         print()
        #     print(ty + " " + str(index) + " " + str(value))


        # To rebuild the model, run operations in reverse order
        operations.reverse()

        # Write the result in output file
        output_model = obja.Output(output, random_color=False)

        for (ty, index, value) in operations:
            if ty == 'v':
                output_model.add_vertex(index, value)
            elif ty == 'f':
                output_model.add_face(index, value)
            elif ty == 'ef':
                output_model.edit_face(index, value)
            elif ty[0] == '#':
                pass
            else:
                raise RuntimeError('Operation non reconnue')


    def compute_K(self):
        '''
        Calcule l'erreur quadratique fondamentale K pour chaque plan associé à chaque face
        '''
        K = {}

        for (face_index, face) in self.batch_faces:
            # Calcul des paramètres du plan de la face : ax + by + cz = d

            p1, p2, p3 = self.vertices[face.a], self.vertices[face.b], self.vertices[face.c]
            e1 = p3 - p1
            e2 = p2 - p1

            cp = np.cross(e1, e2)
            a, b, c = cp

            d = np.dot(cp, p3)
            K[face_index] = np.array([[a*a, a*b, a*c, a*d],
                                      [a*b, b*b, b*c, b*d],
                                      [a*c, b*c, c*c, c*d],
                                      [a*d, b*d, c*d, d*d]])

        return K


    def compute_error(self):    
        '''
        Calcule l'erreur induite sur le modele pour chaque arete retiree
        '''
        edges_error = {} # dic qui à la paire (v1,v2) associe l'erreur, avec ind_v1 < ind_v2

        for (ind_v1, v1) in self.batch_vertices:
            for (face_index, face) in self.batch_faces:
                if ind_v1 in [face.a,face.b,face.c]:
                    # if other indices of the face are greater than ind_v1, we consider the edge
                    for ind_v2 in [face.a,face.b,face.c]:
                        if ind_v2 > ind_v1:
                            # calculer l'erreur de (v1,v2) lorsqu'on retire v1
                            v1_homogene = np.concatenate((v1, [1]),0)
                            erreur = np.matmul(v1_homogene, np.matmul((self.vertex_error(ind_v1) + self.vertex_error(ind_v2)), v1_homogene))
                            
                            edges_error[(ind_v1,ind_v2)] = erreur

        return edges_error 


    def vertex_error(self, ind_v):
        '''
        Calcule la matrice Q qui correspond à la somme des K des surfaces voisines
        '''
        Q = np.zeros((4,4))
        nb_faces = 0

        for (face_index, face) in self.batch_faces:
            if ind_v in [face.a,face.b,face.c]:
                Q += self.K[face_index]
                nb_faces += 1

        return Q / nb_faces


    def select_valid_edges(self, edges, n):
        '''
        Selectionne les n edges ayant le moins d'impact et correspondant aux trois conditions topologiques de suppression
        '''
        if len(self.batch_vertices) < n:
            n = len(self.batch_vertices)

        removed_vertices = set()
        removed_edges = set()

        nb_valid_edges = 0

        for e in edges:
            # Si on a atteint le nombre d'aretes necessaires, on sort
            if nb_valid_edges == n:
                break

            (v1,v2) = e

            v1_neighbors = set()
            v2_neighbors = set()
            shared_neighbor_faces = set()

            condition1 = not (v1 in removed_vertices or v2 in removed_vertices)

            if condition1:       
                # On recupere les voisins des sommets
                for batch_index, (face_index,face) in enumerate(self.batch_faces):
                    face_vertices = [face.a,face.b,face.c]

                    v1_in_face = False
                    if v1 in face_vertices:
                        v1_neighbors.update(face_vertices)
                        v1_in_face = True

                    if v2 in face_vertices:
                        v2_neighbors.update(face_vertices)
                        if v1_in_face:
                            shared_neighbor_faces.add(batch_index)

                condition3 = True
                for (w1,w2) in removed_edges:
                    if (w1 in v1_neighbors and w2 in v2_neighbors) or (w1 in v2_neighbors and w2 in v1_neighbors):
                        condition3 = False
                        break

                if condition3:  
                    shared_neighbors = v1_neighbors.intersection(v2_neighbors)

                    condition2 = True
                    for w in shared_neighbors:
                        is_w_inside_shared_face = False
                        # On regarde si w appartient à l'une des faces communes à v1 et v2
                        for batch_index in shared_neighbor_faces:
                            face = self.batch_faces[batch_index][1]

                            # Si c'est le cas, on sort
                            if w in [face.a,face.b,face.c]:
                                is_w_inside_shared_face = True
                                break
                        
                        # Si ce n'est pas le cas, la condition n'est pas respectée
                        if not is_w_inside_shared_face:
                            condition2 = False
                            break

                    if condition2:                              
                        # Si les trois conditions sont reunies
                        removed_vertices.update([v1,v2])
                        removed_edges.add(e)
                        nb_valid_edges += 1
            
        return removed_edges



def main():
    """
    Runs the program on the model given as parameter.
    """
    
    if len(sys.argv) != 2:
        filename = input('Obj file name? (without extension)')
    else:
        filename = sys.argv[1]
        
    model_path = 'models/'
    output_path = 'output/'

    np.seterr(invalid = 'raise')
    model = Decimater()
    model.parse_file(model_path + filename + '.obj')

    ind = 1
    while isfile(output_path + filename + '_squeeze' + str(ind) + '.obja'):
        ind += 1

    with open(output_path + filename + '_squeeze' + str(ind) + '.obja', 'w') as output_fs:
        model.contract(output_fs)


if __name__ == '__main__':
    main()
