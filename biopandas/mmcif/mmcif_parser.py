import gzip
import re

try:
    import json
except:
    import simplejson as json

try:
    str.partition

    def partitionString(string, sep):
        return string.partition(sep)

except:

    def partitionString(string, sep):
        tmp = string.split(sep)
        return [tmp.pop(0), sep, sep.join(tmp)]


class _loop:
    def __init__(self, parserObj):
        self.parserObj = parserObj
        self.length = 0
        self.refID = -1

        self.refList = []

        self.namesDefined = False

    def addName(self, name):
        catName = type(name) == str and partitionString(name, ".") or ["", "", ""]
        if catName[1]:
            if not catName[0] in self.parserObj.currentTarget[-2]:
                self.parserObj.currentTarget[-2][catName[0]] = {}
            if not catName[2] in self.parserObj.currentTarget[-2][catName[0]]:
                self.parserObj.currentTarget[-2][catName[0]][catName[2]] = []
            self.refList.append(
                self.parserObj.currentTarget[-2][catName[0]][catName[2]]
            )
        else:
            if not catName[0] in self.parserObj.currentTarget[-2]:
                self.parserObj.currentTarget[-2][catName[0]] = []
            self.refList.append(self.parserObj.currentTarget[-2][catName[0]])
        self.length = len(self.refList)

    def pushValue(self, value):
        if not self.namesDefined:
            self.namesDefined = True
        target = self.nextTarget()
        if value == "stop_":
            return self.stopPush()
        target.append(value)

    def nextTarget(self):
        self.refID = (self.refID + 1) % self.length
        return self.refList[self.refID]

    def stopPush(self):
        self.refID = -1


def specialSplit(content):
    output = [["", False]]
    quote = False
    length = len(content)
    log = ""
    for c in range(length):
        isWS = content[c] in [" ", "\t"]
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
        elif not quote and isWS and output[-1][0] != "":
            output.append(["", False])
        elif not quote and content[c] == "#":
            break
        elif not isWS or quote:
            output[-1][0] += content[c]
            output[-1][1] = quote
    if output[-1][0] == "":
        output.pop()
    return output


class targetSetter:
    def __init__(self, obj, key):
        self.obj = obj
        self.key = key

    def setValue(self, value):
        self.obj[self.key] = value


class CIFparser:
    def __init__(self):
        self.data = {}
        self.currentTarget = None
        self.loopPointer = None

    def parseString(self, contents):
        multi_line_mode = False
        buffer = []
        for line in contents.splitlines():
            Z = line[:1]
            line = line.strip()
            if Z == ";":
                if multi_line_mode:
                    self.setDataValue("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.processContent(specialSplit(line))

    def parse(self, fileobj):
        multi_line_mode = False
        buffer = []
        for line in fileobj.readlines():
            Z = line[:1]
            line = line.strip()
            if Z == ";":
                if multi_line_mode:
                    self.setDataValue("\n".join(buffer))
                else:
                    buffer = []
                multi_line_mode = not multi_line_mode
                line = line[1:].strip()
            if multi_line_mode:
                buffer.append(line)
            else:
                self.processContent(specialSplit(line))

    def processContent(self, content):
        for c, quoted in content:
            if c == "global_" and not quoted:
                self.loopPointer = None
                self.selectGlobal()
            elif c[:5] == "data_" and not quoted:
                self.loopPointer = None
                self.selectData(c)
            elif c[:5] == "save_" and not quoted:
                self.loopPointer = None
                if c[5:]:
                    self.selectFrame(c)
                else:
                    self.endFrame()
            elif c == "loop_" and not quoted:
                self.loopPointer = _loop(self)
            elif c[:1] == "_" and not quoted:
                self.setDataName(c[1:])
            else:
                self.setDataValue(c)

    def setDataName(self, name):
        if self.loopPointer != None:
            if self.loopPointer.namesDefined:
                self.loopPointer = None
            else:
                return self.loopPointer.addName(name)
        name = partitionString(name, ".")
        self.currentTarget.pop()
        if name[1]:
            if name[0] not in self.currentTarget[-1]:
                self.currentTarget[-1][name[0]] = {}
            self.currentTarget[-1][name[0]][name[2]] = ""
            self.currentTarget = self.currentTarget + [
                targetSetter(self.currentTarget[-1][name[0]], name[2])
            ]
        else:
            self.currentTarget[-1][name[0]] = ""
            self.currentTarget = self.currentTarget + [
                targetSetter(self.currentTarget[-1], name[0])
            ]

    def setDataValue(self, value):
        if self.loopPointer is None:
            self.currentTarget[-1].setValue([value])

        else:
            self.loopPointer.pushValue(value)

    def selectGlobal(self):
        self.currentTarget = [self.data, self.data, None]

    def selectData(self, name):
        if name not in self.data:
            self.data[name] = {}
        self.currentTarget = [self.data, self.data[name], None]

    def selectFrame(self, name=""):
        if name not in self.currentTarget[1]:
            self.currentTarget[1][name] = {}
        self.currentTarget = self.currentTarget[:2] + [
            self.currentTarget[1][name],
            None,
        ]

    def endData(self):
        self.currentTarget = self.currentTarget[:2]

    def endFrame(self):
        self.currentTarget = self.currentTarget[:3]


####################################################################################################################################################


class __CIFfloat__(float):
    def __repr__(self):
        return "%.15g" % self


class __CIFint__(int):
    def __repr__(self):
        return str(self)


def __CIFfloatRange__(inp):
    try:
        pos = inp.index("-", 1)
        return (__CIFfloat__(inp[:pos]), __CIFfloat__(inp[pos + 1 :]))
    except:
        return (__CIFfloat__(inp),)


def __CIFintRange__(inp):
    try:
        pos = inp.index("-", 1)
        return (__CIFint__(inp[:pos]), __CIFint__(inp[pos + 1 :]))
    except:
        return (__CIFint__(inp),)


def __loadCIFdic__(dicFile, force=False):
    jsfDic = f"{dicFile[:-4]}.json"
    jsf = f"{dicFile[:-4]}_summary.json"
    dic = {}
    try:
        if force:
            throw
        dic = json.loads(open(jsf).read())
    except:
        parser = CIFparser()
        parser.parse(open(dicFile))
        json.dump(parser.data, open(jsfDic, "w"))
        for k, v in parser.data["data_mmcif_pdbx.dic"].items():
            if type(v) != dict or "item_type" not in v:
                continue
            name = partitionString(k[6:], ".")
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
                typing[k][k2] = __CIFfloat__
            elif v2 == "float-range":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __CIFfloatRange__
            elif v2 == "int":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __CIFint__
            elif v2 == "int-range":
                if k not in typing:
                    typing[k] = {}
                typing[k][k2] = __CIFintRange__
    return typing


def __dumpCIF__(jso):
    return __dumpPart__(jso)


__cifStrCheck__ = re.compile(r"[\\s\(\)]")
__cifStrNLCheck__ = re.compile(r"[\n]")


def __dumpStr__(inp):
    if inp == None:
        return "?"
    else:
        if type(inp) != str and type(inp) != unicode:
            return str(inp)
        if re.search(__cifStrNLCheck__, inp) != None:
            return "\n;%s\n;" % inp
        if re.search(__cifStrCheck__, inp) != None:
            return "'%s'" % inp
        else:
            return inp


def __padString__(inp, flength):
    return inp + (" " * (flength - len(inp)))


def __dumpCat__(k, v):
    output = "#\n"
    noi = len(v[v.keys()[0]])
    if noi == 1:
        pad = 0
        for k2 in v.keys():
            if len(k2) > pad:
                pad = len(k2)
        pad += 3
        for k2 in v.keys():
            output += "_%s.%s%s\n" % (k, __padString__(k2, pad), __dumpStr__(v[k2][0]))
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
            for k2 in v.keys():
                tmp2.append(__dumpStr__(v[k2][i]))

        for j in range(len(tmp1[0])):
            pad = 0
            for i in range(len(tmp1)):
                if tmp1[i][j][:2] != "\n;" and len(tmp1[i][j]) > pad:
                    pad = len(tmp1[i][j])
            pad += 1
            for i in range(len(tmp1)):
                if tmp1[i][0][:2] != "\n;":
                    tmp1[i][j] = __padString__(tmp1[i][j], pad)

        for i in range(noi):
            output += "".join(tmp1[i]) + "\n"

    return output.strip() + "\n"


def __dumpPart__(jso):
    inner = True
    output = ""
    for k, v in jso.items():
        if isinstance(v, dict):
            if k[:5] != "data_" and k[:5] != "save_" and k[:7] != "global_":
                output += __dumpCat__(k, v)
            else:
                output += k + "\n"
                output += __dumpPart__(v)
                inner = False
    if inner:
        return output + "#\n"
    else:
        return output


def __loadCIFData__(data, doClean=True, doType=True):
    parser = CIFparser()
    if type(data) == str:
        parser.parseString(data)
    else:
        parser.parse(data)  # fileobj

    if not doClean:
        return parser.data

    for k, v in parser.data.items():
        for k2, v2 in v.items():
            for k3, v3 in v2.items():
                for i in range(len(v3)):
                    v2[k3][i] = v3[i] not in ["?", "."] and v3[i] or None

    if not doType or not __mmcifTyping__:
        return parser.data

    for struct, data in parser.data.items():
        for k, v in __mmcifTyping__.items():
            if k not in data:
                continue
            for k2, v2 in v.items():
                if k2 in data[k]:
                    for r in range(len(data[k][k2])):
                        try:
                            data[k][k2][r] = v2(data[k][k2][r])
                        except:
                            pass

    return parser.data


def __loadCIF__(cifFile, doClean=True, doType=True):
    parser = CIFparser()
    if cifFile[-3:].lower() == ".gz":
        parser.parse(gzip.open(cifFile))
    else:
        parser.parse(open(cifFile))

    if not doClean:
        return parser.data

    for k, v in parser.data.items():
        for k2, v2 in v.items():
            for k3, v3 in v2.items():
                for i in range(len(v3)):
                    v2[k3][i] = v3[i] not in ["?", "."] and v3[i] or None

    if not doType or not __mmcifTyping__:
        return parser.data

    for struct, data in parser.data.items():
        for k, v in __mmcifTyping__.items():
            if k not in data:
                continue
            for k2, v2 in v.items():
                if k2 in data[k]:
                    for r in range(len(data[k][k2])):
                        try:
                            data[k][k2][r] = v2(data[k][k2][r])
                        except:
                            pass

    return parser.data


__mmcifTyping__ = None
