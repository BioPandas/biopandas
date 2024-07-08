# the cif parser relieves heavily on the parser from pdb_japan
# (author: gert-jan bekker).
# python cif parser: https://gitlab.com/pdbjapan/tools/cif-parsers
# license: mit
# see https://gitlab.com/pdbjapan/tools/cif-parsers/blob/master/license
import contextlib
import gzip
import json
import re


def partition_string(string, sep):
    return string.partition(sep)


class LoopMMCIF:
    def __init__(self, parser_obj):
        self.parser_obj = parser_obj
        self.length = 0
        self.ref_id = -1
        self.ref_list = []
        self.names_defined = False

    def add_name(self, name):
        cat_name = (
            isinstance(name, str) and partition_string(name, ".") or ["", "", ""]
        )
        if cat_name[1]:
            if cat_name[0] not in self.parser_obj.current_target[-2]:
                self.parser_obj.current_target[-2][cat_name[0]] = {}
            if (
                cat_name[2]
                not in self.parser_obj.current_target[-2][cat_name[0]]
            ):
                self.parser_obj.current_target[-2][cat_name[0]][
                    cat_name[2]
                ] = []
            self.ref_list.append(
                self.parser_obj.current_target[-2][cat_name[0]][cat_name[2]]
            )
        else:
            if cat_name[0] not in self.parser_obj.current_target[-2]:
                self.parser_obj.current_target[-2][cat_name[0]] = []
            self.ref_list.append(
                self.parser_obj.current_target[-2][cat_name[0]]
            )
        self.length = len(self.ref_list)

    def push_value(self, value):
        if not self.names_defined:
            self.names_defined = True
        target = self.next_target()
        if value == "stop_":
            return self.stop_push()
        target.append(value)

    def next_target(self):
        self.ref_id = (self.ref_id + 1) % self.length
        return self.ref_list[self.ref_id]

    def stop_push(self):
        self.ref_id = -1


def special_split(content):
    output = [["", False]]
    quote = False
    length = len(content)
    for c in range(length):
        is_ws = content[c] in [" ", "\t"]
        if content[c] in ["'", '"'] and (
            (
                c == 0
                or content[c - 1] == " "
                or content[c - 1] == "\t"
                or c == length - 1
                or content[c + 1] == " "
                or content[c + 1] == "\t"
            )
        ):
            quote = not quote
        elif not quote and is_ws and output[-1][0] != "":
            output.append(["", False])
        elif not quote and content[c] == "#":
            break
        elif not is_ws or quote:
            output[-1][0] += content[c]
            output[-1][1] = quote
    if output[-1][0] == "":
        output.pop()
    return output


class TargetSetter:
    def __init__(self, obj, key):
        self.obj = obj
        self.key = key

    def set_value(self, value):
        self.obj[self.key] = value


class CIFParser:
    def __init__(self):
        self.data = {}
        self.current_target = None
        self.loop_pointer = None

    def parse_string(self, contents):
        multi_line_mode = False
        buffer = []
        for line in contents.splitlines():
            z = line[:1]
            line = line.strip()
            if z == ";":
                if multi_line_mode:
                    self.set_data_value("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.process_content(special_split(line))

    def parse(self, fileobj):
        multi_line_mode = False
        buffer = []
        for line in fileobj.readlines():
            z = line[:1]
            line = line.strip()
            if z == ";":
                if multi_line_mode:
                    self.set_data_value("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.process_content(special_split(line))

    def process_content(self, content):
        for c, quoted in content:
            if c == "global_" and not quoted:
                self.loop_pointer = None
                self.select_global()
            elif c[:5] == "data_" and not quoted:
                self.loop_pointer = None
                self.select_data(c)
            elif c[:5] == "save_" and not quoted:
                self.loop_pointer = None
                if c[5:]:
                    self.select_frame(c)
                else:
                    self.end_frame()
            elif c == "loop_" and not quoted:
                self.loop_pointer = LoopMMCIF(self)
            elif c[:1] == "_" and not quoted:
                self.set_data_name(c[1:])
            else:
                self.set_data_value(c)

    def set_data_name(self, name):
        if self.loop_pointer is not None:
            if self.loop_pointer.names_defined:
                self.loop_pointer = None
            else:
                return self.loop_pointer.add_name(name)
        name = partition_string(name, ".")
        self.current_target.pop()
        if name[1]:
            if name[0] not in self.current_target[-1]:
                self.current_target[-1][name[0]] = {}
            self.current_target[-1][name[0]][name[2]] = ""
            self.current_target = self.current_target + [
                TargetSetter(self.current_target[-1][name[0]], name[2])
            ]
        else:
            self.current_target[-1][name[0]] = ""
            self.current_target = self.current_target + [
                TargetSetter(self.current_target[-1], name[0])
            ]

    def set_data_value(self, value):
        if self.loop_pointer is None:
            self.current_target[-1].set_value([value])

        else:
            self.loop_pointer.push_value(value)

    def select_global(self):
        self.current_target = [self.data, self.data, None]

    def select_data(self, name):
        if name not in self.data:
            self.data[name] = {}
        self.current_target = [self.data, self.data[name], None]

    def select_frame(self, name=""):
        if name not in self.current_target[1]:
            self.current_target[1][name] = {}
        self.current_target = self.current_target[:2] + [
            self.current_target[1][name],
            None,
        ]

    def end_data(self):
        self.current_target = self.current_target[:2]

    def end_frame(self):
        self.current_target = self.current_target[:3]


class __CIFFloat__(float):
    def __repr__(self):
        return "%.15g" % self


class __CIFInt__(int):
    def __repr__(self):
        return str(self)


def __cif_float_range__(inp):
    try:
        pos = inp.index("-", 1)
        return (__CIFFloat__(inp[:pos]), __CIFFloat__(inp[pos + 1:]))
    except Exception:
        return (__CIFFloat__(inp),)


def __cif_int_range__(inp):
    try:
        pos = inp.index("-", 1)
        return (__CIFInt__(inp[:pos]), __CIFInt__(inp[pos + 1:]))
    except Exception:
        return (__CIFInt__(inp),)


def __load_cif_dic__(dic_file, force=False):
    jsf_dic = f"{dic_file[:-4]}.json"
    jsf = f"{dic_file[:-4]}_summary.json"
    dic = {}
    try:
        if force:
            throw
        dic = json.loads(open(jsf).read())
    except Exception:
        parser = CIFParser()
        parser.parse(open(dic_file))
        json.dump(parser.data, open(jsf_dic, "w"))
        for k, v in parser.data["data_mmcif_pdbx.dic"].items():
            if not isinstance(v, dict) or "item_type" not in v:
                continue
            name = partition_string(k[6:], ".")
            if name[0] not in dic:
                dic[name[0]] = {}
            dic[name[0]][name[2]] = v["item_type"]["code"][0].strip()
        json.dump(dic, open(jsf, "w"))

    typing = {}
    for k, v in dic.items():
        for k2, v2 in v.items():
            if v2 == "float":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __CIFFloat__
            elif v2 == "float-range":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __cif_float_range__
            elif v2 == "int":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __CIFInt__
            elif v2 == "int-range":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __cif_int_range__
    return typing


def __dump_cif__(jso):
    return __dump_part__(jso)


__CIF_STR_CHECK__ = re.compile(r"[\\s\(\)]")
__CIF_STR_NL_CHECK__ = re.compile(r"[\n]")


def __dump_str__(inp):
    if inp is None:
        return "?"
    if not isinstance(inp, str):
        return str(inp)
    if re.search(__CIF_STR_NL_CHECK__, inp) is not None:
        return "\n;%s\n;" % inp
    return (
        "'%s'" % inp if re.search(__CIF_STR_CHECK__, inp) is not None else inp
    )


def __pad_string__(inp, flength):
    return inp + (" " * (flength - len(inp)))


def __dump_cat__(k, v):
    output = "#\n"
    noi = len(v[v.keys()[0]])
    if noi == 1:
        pad = 0
        for k2 in v.keys():
            if len(k2) > pad:
                pad = len(k2)
        pad += 3
        for k2 in v.keys():
            output += "_%s.%s%s\n" % (
                k,
                __pad_string__(k2, pad),
                __dump_str__(v[k2][0]),
            )
    else:
        output += "loop_\n"
        pad = []
        for k2 in v.keys():
            output += "_%s.%s\n" % (k, k2)
            pad.append(0)
        tmp1 = []
        for i in range(noi):
            tmp2 = []
            tmp1.append(tmp2)
            tmp2.extend(__dump_str__(v[k2][i]) for k2 in v.keys())
        for j in range(len(tmp1[0])):
            pad = 0
            for item_ in tmp1:
                if item_[j][:2] != "\n;" and len(item_[j]) > pad:
                    pad = len(item_[j])
            pad += 1
            for item in tmp1:
                if item[0][:2] != "\n;":
                    item[j] = __pad_string__(item[j], pad)

        for i in range(noi):
            output += "".join(tmp1[i]) + "\n"

    return output.strip() + "\n"


def __dump_part__(jso):
    inner = True
    output = ""
    for k, v in jso.items():
        if isinstance(v, dict):
            if k[:5] != "data_" and k[:5] != "save_" and k[:7] != "global_":
                output += __dump_cat__(k, v)
            else:
                output += k + "\n"
                output += __dump_part__(v)
                inner = False
    return output + "#\n" if inner else output


def load_cif_data(data, do_clean=True, do_type=True):
    parser = CIFParser()
    if isinstance(data, str):
        parser.parse_string(data)
    else:
        parser.parse(data)  # fileobj

    if not do_clean:
        return parser.data

    for k, v in parser.data.items():
        for k2, v2 in v.items():
            for k3, v3 in v2.items():
                for i in range(len(v3)):
                    v2[k3][i] = v3[i] not in ["?", "."] and v3[i] or None

    if not do_type or not __MMCIF_TYPING__:
        return parser.data

    for _, data in parser.data.items():
        for k, v in __MMCIF_TYPING__.items():
            if k not in data:
                continue
            for k2, v2 in v.items():
                if k2 in data[k]:
                    for r in range(len(data[k][k2])):
                        with contextlib.suppress(KeyError):
                            data[k][k2][r] = v2(data[k][k2][r])
    return parser.data


def __load_cif__(cif_file, do_clean=True, do_type=True):
    parser = CIFParser()
    if cif_file[-3:].lower() == ".gz":
        parser.parse(gzip.open(cif_file))
    else:
        parser.parse(open(cif_file))

    if not do_clean:
        return parser.data

    for k, v in parser.data.items():
        for k2, v2 in v.items():
            for k3, v3 in v2.items():
                for i in range(len(v3)):
                    v2[k3][i] = v3[i] not in ["?", "."] and v3[i]  # or None

    if not do_type or not __MMCIF_TYPING__:
        return parser.data

    for _, data in parser.data.items():
        for k, v in __MMCIF_TYPING__.items():
            if k not in data:
                continue
            for k2, v2 in v.items():
                if k2 in data[k]:
                    for r in range(len(data[k][k2])):
                        with contextlib.suppress(KeyError):
                            data[k][k2][r] = v2(data[k][k2][r])
    return parser.data


__MMCIF_TYPING__ = None
