from rdkit.Chem.Draw import rdMolDraw2D
import tkinter
from io import BytesIO
from typing import Union, Tuple, Optional, List
from PIL import Image, ImageTk
from rdkit import Chem
import threading
import copy
import time


class MolDraw:
    def __init__(
            self,
            size: Tuple[int, int] = (500, 500),
            title: str = 'RDKit Molecule',
            rotate: float = 0,
            show_atom_index: bool = False,
            show_bond_index: bool = False,
            show_stereo_label: bool = False,
            close_time: float = -1,
            stay_in_front: bool = True,
            high_light_atoms: Optional[List[int]] = None,
            high_light_bonds: Optional[List[int]] = None,
            use_plt: bool = False,
            clear_map: bool = True,
            cairo: bool = False,
            use_acs: bool = False,
            save_path:Union[None,str] = None
    ):
        """
        Initialize molecule drawer with specified options.

        Args:
            size: Dimensions of the output image (width, height)
            title: Title for the display window
            rotate: Rotation angle in degrees
            show_atom_index: Whether to display atom indices
            show_bond_index: Whether to display bond indices
            show_stereo_label: Whether to show stereo chemistry labels
            close_time: Time in seconds before window auto-closes (-1 for manual close)
            stay_in_front: Keep window on top of others
            high_light_atoms: List of atom indices to highlight
            high_light_bonds: List of bond indices to highlight
            use_plt: Use matplotlib for display (else Tkinter)
            clear_map: Clear atom maps before drawing
            cairo: Use Cairo backend (else SVG)
            use_acs: Use ACS 1996 drawing style
        """
        self.figure_size = size
        self.title = title
        self.rotate = rotate
        self.show_atom_index = show_atom_index
        self.show_bond_index = show_bond_index
        self.show_stereo_label = show_stereo_label
        self.close_time = close_time
        self.stay_in_front = stay_in_front
        self.high_light_atoms = high_light_atoms or []  # 空和None都是false 从左往右检测 返回[]
        self.high_light_bonds = high_light_bonds or []
        self.use_plt = use_plt
        self.clear_map = clear_map
        self.cairo = cairo
        self.acs = use_acs
        self.d2d = None
        self.save_path = save_path

    def mol_draw(self, mol):
        """Draw RDKit molecule with configured options.

        Args:
            mol: RDKit molecule to draw

        Returns:
            Depending on display method, returns drawer and/or display objects
        """
        if mol is None:
            raise ValueError("Molecule cannot be None")

        plot_mol = copy.deepcopy(mol)
        Chem.RemoveHs(plot_mol)

        if self.cairo:
            d2d = rdMolDraw2D.MolDraw2DCairo(*self.figure_size)
        else:
            d2d = rdMolDraw2D.MolDraw2DSVG(*self.figure_size)

        self._configure_draw_options(d2d.drawOptions())

        d2d.DrawMolecule(
            plot_mol,
            highlightAtoms=self.high_light_atoms,
            highlightBonds=self.high_light_bonds
        )
        d2d.FinishDrawing()
        self.d2d = d2d
        if not self.use_plt:
            self._show_tkinter()
        else:
            self._show_matplotlib()
        if self.save_path:
            self._save_plot()


    def _configure_draw_options(self, opts) -> None:
        """Configure drawing options based on instance settings."""
        if self.acs:
            opts.SetACS1996Mode()
        else:
            opts.addAtomIndices = self.show_atom_index
            opts.addBondIndices = self.show_bond_index
            opts.addStereoAnnotation = self.show_stereo_label
            opts.rotate = self.rotate
            opts.dummyIsotopeLabels = True

    def _show_tkinter(self) -> tkinter.Tk:
        """Display molecule using Tkinter."""
        try:
            sio = BytesIO(self.d2d.GetDrawingText())
            img = Image.open(sio)

            tk_root = tkinter.Tk()
            tk_root.title(self.title)
            tk_pi = ImageTk.PhotoImage(img)

            tk_label = tkinter.Label(tk_root, image=tk_pi)
            tk_label.place(x=0, y=0, width=img.size[0], height=img.size[1])

            tk_root.geometry(f'{img.size[0]}x{img.size[1]}')
            tk_root.lift()

            if self.stay_in_front:
                tk_root.attributes('-topmost', True)

            if self.close_time > 0:
                threading.Thread(
                    target=self._close_window_after_delay,
                    args=(tk_root, self.close_time),
                    daemon=True
                ).start()

            tk_root.mainloop()
            return tk_root

        except Exception as e:
            raise RuntimeError(f"Failed to display molecule with Tkinter: {str(e)}")

    @staticmethod
    def _close_window_after_delay(window, delay: float) -> None:
        """Close window after specified delay."""
        time.sleep(delay)
        window.quit()
        window.destroy()

    def _show_matplotlib(self) -> None:
        """Display molecule using matplotlib."""
        raise NotImplementedError("Matplotlib display not yet implemented")

    def _save_plot(self):
        if self.cairo:
            with open(f'{self.save_path}.png', 'wb') as f:
                f.write(self.d2d.GetDrawingText())
        else:
            with open(f'{self.save_path}.svg', 'w') as f:
                f.write(self.d2d.GetDrawingText())



# def show_mol_with_matplotlib(mol: Union[Mol, None], size: Tuple[int, int] = (200, 200), title: str = 'RDKit Molecule',
#                              svg: bool = False, rotate: float = 0, show_atom: bool = False,
#                              show_bond: bool = False, high_light_atoms: list = [],
#                              stay_in_front=True, show_stereo: bool = False, auto_close: bool = True,
#                              close_time: int = 2, useTemp: Union[bool, str] = False, clearMap: str = True,
#                              *args, **kwargs):
#     """
#         Simply show moleculmes
#     :param mol: rdkit molecule object
#     :param auto_close:
#     :param savePath:
#     :param clearMap: clear map num default True
#     :param show_bond: show bond idx num
#     :param useTemp: use a template to show
#     :param close_time: auto close ticker window time unit second
#     :param rotate: pic rotate degree
#     :param show_stereo: show atom CIP(R/S) in pic
#     :param stay_in_front: window show in front
#     :param show_atom: show atom index in pic
#     :param title: pic title
#     :param size: pic size (ixx,ixy)
#     :return:
#     """
#     plotMol = copy.deepcopy(mol)  # 使用mapnum绘图时有可能影响序号 所以深拷贝一个出来
#     if type(useTemp) is str:
#         try:
#             tempMol = Chem.MolFromSmarts(useTemp)
#             aChem.Compute2DCoords(tempMol)
#             aChem.GenerateDepictionMatching2DStructure(plotMol, tempMol)
#         except Exception as E:
#             print('Plot with temp failed')
#     else:
#         pass
#     if clearMap:
#         for atom in plotMol.GetAtoms():
#             atom.SetAtomMapNum(0)
#
#     Chem.RemoveHs(plotMol)
#     if svg:
#         d2d = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
#     else:
#         d2d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
#     d2d.drawOptions().addAtomIndices = show_atom
#     d2d.drawOptions().addBondIndices = show_bond
#     d2d.drawOptions().addStereoAnnotation = show_stereo
#     d2d.drawOptions().rotate = rotate
#     d2d.drawOptions().bondLineWidth = 0.5
#     d2d.drawOptions().dummyIsotopeLabels = True
#     d2d.DrawMolecule(plotMol, highlightAtoms=high_light_atoms)
#     d2d.FinishDrawing()
#     if not svg:
#         sio = BytesIO(d2d.GetDrawingText())
#         img = Image.open(sio)
#
#         # 使用 matplotlib 显示图像
#         plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)  # 调整尺寸（像素转英寸）
#         plt.imshow(img)
#         plt.axis('off')  # 关闭坐标轴
#         plt.title(title)
#         plt.show()
#     else:
#         with open('temp.svg', 'w') as f:
#             f.write(d2d.GetDrawingText())
#
#     return d2d
#
#
# def show_reaction(reaction: Union[Mol, None], size: Tuple[int, int] = (1000, 1000), title: str = 'RDKit Molecule',
#                   rotateDegree: float = 0, showAtom: bool = False, stayInFront=True, showStereoAnnotation: bool = False,
#                   autoclose: int = 1.5, *args, **kwargs):
#     """
#         Simply show molecules
#     :param autoclose:
#     :param reaction:
#     :param rotateDegree: pic rotare degree
#     :param showStereoAnnotation: show atom CIP(R/S) in pic
#     :param stayInFront:window show in front
#     :param showAtom: show atom index in pic
#     :param title: pic title
#     :param mol: rdkit molecule object
#     :param size: pic size (ixx,ixy)
#     :return:
#     """
#     d2d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
#     d2d.drawOptions().addAtomIndices = False
#     d2d.DrawReaction(reaction, **kwargs)
#     d2d.FinishDrawing()
#     sio = BytesIO(d2d.GetDrawingText())
#     img = Image.open(sio)
#     tkRoot = tkinter.Tk()
#     tkRoot.title(title)
#     tkPI = ImageTk.PhotoImage(img)
#     tkLabel = tkinter.Label(tkRoot, image=tkPI)
#     tkLabel.place(x=0, y=0, width=img.size[0], height=img.size[1])
#     tkRoot.geometry('%dx%d' % img.size)
#     tkRoot.lift()
#     if stayInFront:
#         tkRoot.attributes('-topmost', True)
#
#     def close_window(auto_close):
#         tkRoot.update()
#         time.sleep(auto_close)
#         tkRoot.quit()
#         tkRoot.destroy()
#
#     if autoclose:
#         threading.Thread(target=close_window(autoclose)).start()
#     tkRoot.mainloop()
#     return d2d
#
#
# def save_mol(savePath: str, mol: Union[Mol, None], size: Tuple[int, int] = (500, 500), title: str = 'RDKit Molecule',
#              rotateDegree: float = 0, showAtom: bool = False, showStereoAnnotation: bool = False,
#              *args, **kwargs):
#     """
#
#     :param title:
#     :param showStereoAnnotation:
#     :param size:
#     :param rotateDegree:
#     :param mol:
#     :param text:
#     :param showBond:
#     :param showAtom:
#     :param savePath:
#     :param template:
#     :param save:
#     :return:
#     """
#     # tempMol = Chem.MolFromSmarts(
#     #     '[#6]12~[#6]~[#6]~[#6]3~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~3~[#6]~1~[#6]~[#6]~[#6]~[#6]~2')
#     # aChem.Compute2DCoords(tempMol)
#     # aChem.GenerateDepictionMatching2DStructure(mol, tempMol)
#     d2d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
#     d2d.drawOptions().addAtomIndices = showAtom
#     d2d.drawOptions().addStereoAnnotation = showStereoAnnotation
#     d2d.drawOptions().rotate = rotateDegree
#     d2d.DrawMolecule(mol)
#     d2d.FinishDrawing()
#     svg_text = d2d.GetDrawingText()
#     with open(savePath, 'wb') as png_file:
#         png_file.write(svg_text)
#     png_file.close()
#     return True
